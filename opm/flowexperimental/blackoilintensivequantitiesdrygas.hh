// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::BlackOilIntensiveQuantities
 */
#ifndef EWOMS_BLACK_OIL_INTENSIVE_QUANTITIES_DRYGAS_HH
#define EWOMS_BLACK_OIL_INTENSIVE_QUANTITIES_DRYGAS_HH
#include <opm/models/blackoil/blackoilproperties.hh>
#include <opm/models/blackoil/blackoilsolventmodules.hh>
#include <opm/models/blackoil/blackoilextbomodules.hh>
#include <opm/models/blackoil/blackoilpolymermodules.hh>
#include <opm/models/blackoil/blackoilfoammodules.hh>
#include <opm/models/blackoil/blackoilbrinemodules.hh>
#include <opm/models/blackoil/blackoilenergymodules.hh>
#include <opm/models/blackoil/blackoildiffusionmodule.hh>
#include <opm/models/blackoil/blackoilmicpmodules.hh>
#include <opm/material/fluidstates/BlackOilFluidState.hpp>
#include <opm/material/common/Valgrind.hpp>
#include <opm/input/eclipse/EclipseState/Grid/FaceDir.hpp>
#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/utility/CopyablePtr.hpp>

#include <dune/common/fmatrix.hh>
#include <opm/material/fluidmatrixinteractions/EclDefaultMaterialExperimental.hpp>
#include <cstring>
#include <utility>

#include <fmt/format.h>

namespace Opm {
/*!
 * \ingroup BlackOilModel
 * \ingroup IntensiveQuantities
 *
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the black-oil model.
 */
template <class TypeTag>
class BlackOilIntensiveQuantitiesDryGas
    : public GetPropType<TypeTag, Properties::DiscIntensiveQuantities>
    , public GetPropType<TypeTag, Properties::FluxModule>::FluxIntensiveQuantities
    , public BlackOilDiffusionIntensiveQuantities<TypeTag, getPropValue<TypeTag, Properties::EnableDiffusion>() >
    , public BlackOilSolventIntensiveQuantities<TypeTag>
    , public BlackOilExtboIntensiveQuantities<TypeTag>
    , public BlackOilPolymerIntensiveQuantities<TypeTag>
    , public BlackOilFoamIntensiveQuantities<TypeTag>
    , public BlackOilBrineIntensiveQuantities<TypeTag>
    , public BlackOilEnergyIntensiveQuantities<TypeTag>
    , public BlackOilMICPIntensiveQuantities<TypeTag>
{
    using ParentType = GetPropType<TypeTag, Properties::DiscIntensiveQuantities>;
    using Implementation = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluxModule = GetPropType<TypeTag, Properties::FluxModule>;

    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    enum { enableSolvent = getPropValue<TypeTag, Properties::EnableSolvent>() };
    enum { enableExtbo = getPropValue<TypeTag, Properties::EnableExtbo>() };
    enum { enablePolymer = getPropValue<TypeTag, Properties::EnablePolymer>() };
    enum { enableFoam = getPropValue<TypeTag, Properties::EnableFoam>() };
    enum { enableBrine = getPropValue<TypeTag, Properties::EnableBrine>() };
    //enum { enableEvaporation = getPropValue<TypeTag, Properties::EnableEvaporation>() };
    enum { enableVapwat = getPropValue<TypeTag, Properties::EnableVapwat>() };
    enum { enableSaltPrecipitation = getPropValue<TypeTag, Properties::EnableSaltPrecipitation>() };
    enum { enableTemperature = getPropValue<TypeTag, Properties::EnableTemperature>() };
    enum { enableEnergy = getPropValue<TypeTag, Properties::EnableEnergy>() };
    enum { enableDiffusion = getPropValue<TypeTag, Properties::EnableDiffusion>() };
    enum { enableMICP = getPropValue<TypeTag, Properties::EnableMICP>() };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };
    enum { waterCompIdx = FluidSystem::waterCompIdx };
    enum { oilCompIdx = FluidSystem::oilCompIdx };
    enum { gasCompIdx = FluidSystem::gasCompIdx };
    enum { waterPhaseIdx = FluidSystem::waterPhaseIdx };
    enum { oilPhaseIdx = FluidSystem::oilPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { dimWorld = GridView::dimensionworld };
    enum { compositionSwitchIdx = Indices::compositionSwitchIdx };

    static constexpr bool compositionSwitchEnabled = Indices::compositionSwitchIdx >= 0;
    static constexpr bool waterEnabled = Indices::waterEnabled;
    static constexpr bool gasEnabled = Indices::gasEnabled;
    static constexpr bool oilEnabled = Indices::oilEnabled;

    using Toolbox = MathToolbox<Evaluation>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using FluxIntensiveQuantities = typename FluxModule::FluxIntensiveQuantities;
    using DiffusionIntensiveQuantities = BlackOilDiffusionIntensiveQuantities<TypeTag, enableDiffusion>;

    using DirectionalMobilityPtr = Opm::Utility::CopyablePtr<DirectionalMobility<TypeTag, Evaluation>>;


public:
    using FluidState = BlackOilFluidState<Evaluation,
                                          FluidSystem,
                                          enableTemperature,
                                          enableEnergy,
                                          compositionSwitchEnabled,
                                          enableVapwat,
                                          enableBrine,
                                          enableSaltPrecipitation,
                                          false,
                                          Indices::numPhases>;
    using ScalarFluidState = BlackOilFluidState<Scalar,
                                                FluidSystem,
                                                enableTemperature,
                                                enableEnergy,
                                                compositionSwitchEnabled,
                                                enableVapwat,
                                                enableBrine,
                                                enableSaltPrecipitation,
                                                false,
                                                Indices::numPhases>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    BlackOilIntensiveQuantitiesDryGas()
    {
        if (compositionSwitchEnabled) {
            fluidState_.setRs(0.0);
            fluidState_.setRv(0.0);
        }
        if (enableVapwat) {
            fluidState_.setRvw(0.0);
        }
    }
    BlackOilIntensiveQuantitiesDryGas(const BlackOilIntensiveQuantitiesDryGas& other) = default;

    BlackOilIntensiveQuantitiesDryGas& operator=(const BlackOilIntensiveQuantitiesDryGas& other) = default;

    /*!
     * \copydoc IntensiveQuantities::update
     */
    void update(const ElementContext& elemCtx, unsigned dofIdx, unsigned timeIdx)
    {
        ParentType::update(elemCtx, dofIdx, timeIdx);

        const auto& problem = elemCtx.problem();
        const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        //const auto& linearizationType = problem.model().linearizer().getLinearizationType();
        unsigned globalSpaceIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
        this->update(problem,priVars,globalSpaceIdx,timeIdx);
    }
    void update(const Problem& problem,
                const PrimaryVariables& priVars,
                unsigned globalSpaceIdx,
                unsigned timeIdx)
    {
        OPM_TIMEBLOCK_LOCAL(UpdateIntensiveQuantitiesGenneral);
        const auto& waterpvt = FluidSystem::waterPvt().template getRealPvt<Opm::WaterPvtApproach::ConstantCompressibilityWater>();
        const auto& gaspvt = FluidSystem::gasPvt().template getRealPvt<Opm::GasPvtApproach::DryGas>();
        const auto& oilpvt = FluidSystem::oilPvt().template getRealPvt<Opm::OilPvtApproach::LiveOil>();;
        this->update(problem,priVars, globalSpaceIdx, timeIdx, waterpvt, gaspvt, oilpvt);
    }

    template<class WaterPvtT, class OilPvtT,class GasPvtT>
    void update(const Problem& problem,
                const PrimaryVariables& priVars,
                unsigned globalSpaceIdx,
                unsigned timeIdx,
                const WaterPvtT& waterpvt,
                const GasPvtT& gaspvt,
                const OilPvtT& oilpvt)
    {
        OPM_TIMEBLOCK_LOCAL(UpdateIntensiveQuantities);
        Evaluation T=298.0;
        std::array<bool, numPhases> saturated;
        for(int i=0; i< numPhases; ++i){
            saturated[i] = true;
        }
        std::array<Evaluation, numPhases> viscosity;

        //const auto& priVars = elemCtx.primaryVars(dofIdx, timeIdx);
        const auto& linearizationType = Opm::LinearizationType();//problem.model().linearizer().getLinearizationType();//NB ??
        // Scalar RvMax = FluidSystem::enableVaporizedOil()
        //     ? problem.maxOilVaporizationFactor(timeIdx, globalSpaceIdx)
        //     : 0.0;
        // Scalar RsMax = FluidSystem::enableDissolvedGas()
        //     ? problem.maxGasDissolutionFactor(timeIdx, globalSpaceIdx)
        //     : 0.0;

        Scalar RsMax = 1e99;//problem.maxGasDissolutionFactor(timeIdx, globalSpaceIdx);
        //     : 0.0;
        asImp_().updateTemperature_(problem, priVars, globalSpaceIdx, timeIdx, linearizationType);

        unsigned pvtRegionIdx = priVars.pvtRegionIndex();
        fluidState_.setPvtRegionIndex(pvtRegionIdx);

        //asImp_().updateSaltConcentration_(elemCtx, dofIdx, timeIdx);

        // extract the water and the gas saturations for convenience
        Evaluation Sw = 0.0;
        if constexpr (waterEnabled) {
            if (priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Sw) {
                Sw = priVars.makeEvaluation(Indices::waterSwitchIdx, timeIdx);
            } else if (priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Disabled){
                // water is enabled but is not a primary variable i.e. one phase case
                Sw = 1.0;
            }
        }
        Evaluation Sg = 0.0;
        if constexpr (gasEnabled) {
            if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Sg) {
                Sg = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
            } else if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rv) {
                Sg = 1.0 - Sw;
            } else if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Disabled) {
                if constexpr (waterEnabled) {
                    Sg = 1.0 - Sw; // two phase water + gas
                } else {
                    // one phase case
                    Sg = 1.0;
                }
            }
        }
        Valgrind::CheckDefined(Sg);
        Valgrind::CheckDefined(Sw);

        Evaluation So = 1.0 - Sw - Sg;

        // deal with solvent
        if constexpr (enableSolvent)
            So -= priVars.makeEvaluation(Indices::solventSaturationIdx, timeIdx);

        if (FluidSystem::phaseIsActive(waterPhaseIdx))
            fluidState_.setSaturation(waterPhaseIdx, Sw);

        if (FluidSystem::phaseIsActive(gasPhaseIdx))
            fluidState_.setSaturation(gasPhaseIdx, Sg);

        if (FluidSystem::phaseIsActive(oilPhaseIdx))
            fluidState_.setSaturation(oilPhaseIdx, So);

        //asImp_().solventPreSatFuncUpdate_(elemCtx, dofIdx, timeIdx);

        // now we compute all phase pressures

        std::array<Evaluation, numPhases> pC = {0, 0, 0};
        //computeRelpermAndPC(mobility_, pC, problem, fluidState_, globalSpaceIdx);
        //const auto& materialParams = problem.materialLawParams(globalSpaceIdx);
        //MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);
        //problem.updateRelperms(mobility_, dirMob_, fluidState_, globalSpaceIdx);
        //MaterialLaw::capillaryPressures(pC, materialParams, fluidState);
        //MaterialLaw::relativePermeabilities(mobility, materialParams, fluidState);
        //const auto& materialParams = problem.materialLawParams(globalSpaceIdx);
        //const auto& materialParams = problem.materialLawParams(0);
        {
            OPM_TIMEBLOCK_LOCAL(RelpermEvaluation);
            const auto& materialParams = problem.materialLawParams(0).template getRealParams<Opm::EclMultiplexerApproach::Default>();
            const auto& gasoilparams =materialParams.gasOilParams().drainageParams().effectiveLawParams().template getRealParams<Opm::SatCurveMultiplexerApproach::PiecewiseLinear>();
            const auto& oilwaterparams =materialParams.oilWaterParams().drainageParams().effectiveLawParams().template getRealParams<Opm::SatCurveMultiplexerApproach::PiecewiseLinear>();
            Scalar Swco = materialParams.Swl();
            Opm::EclDefaultMaterialExperimental<typename MaterialLaw::DefaultMaterial>::relativePermeabilities(mobility_, Swco, oilwaterparams, gasoilparams, fluidState_);
            //MaterialLaw::DefaultMaterial::relativePermeabilitiesNew(mobility_, materialParams, fluidState_);
            //MaterialLaw::DefaultMaterial::capillaryPressures(pC, materialParams, fluidState_);
            //makMaterialLaw::DefaultMaterial::relativePermeabilitiesNew(mobility_, materialParams, fluidState_);
            //MaterialLaw::DefaultMaterial::relativePermeabilitiesSimple(mobility_, materialParams, fluidState_);
            // for now
            //MaterialLaw::DefaultMaterial::capillaryPressures(pC, materialParams, fluidState_);
            //problem.updateRelperms(mobility_, dirMob_, fluidState_, globalSpaceIdx);

        }
        //mobility_[waterPhaseIdx] = Sw;
        //mobility_[oilPhaseIdx] = So;
        //mobility_[gasPhaseIdx] = Sg;

        // oil is the reference phase for pressure
        if (priVars.primaryVarsMeaningPressure() == PrimaryVariables::PressureMeaning::Pg) {
            const Evaluation& pg = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                if (FluidSystem::phaseIsActive(phaseIdx))
                    fluidState_.setPressure(phaseIdx, pg + (pC[phaseIdx] - pC[gasPhaseIdx]));
        } else if (priVars.primaryVarsMeaningPressure() == PrimaryVariables::PressureMeaning::Pw) {
            const Evaluation& pw = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx);
           for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                if (FluidSystem::phaseIsActive(phaseIdx))
                    fluidState_.setPressure(phaseIdx, pw + (pC[phaseIdx] - pC[waterPhaseIdx]));
        } else {
            //assert(FluidSystem::phaseIsActive(oilPhaseIdx));
            const Evaluation& po = priVars.makeEvaluation(Indices::pressureSwitchIdx, timeIdx);
            for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
                if (FluidSystem::phaseIsActive(phaseIdx))
                    fluidState_.setPressure(phaseIdx, po + (pC[phaseIdx] - pC[oilPhaseIdx]));
        }

        // update the Saturation functions for the blackoil solvent module.
        //asImp_().solventPostSatFuncUpdate_(elemCtx, dofIdx, timeIdx);

        // update extBO parameters
        //asImp_().zFractionUpdate_(elemCtx, dofIdx, timeIdx);

        Evaluation SoMax = 0.0;
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            SoMax = max(fluidState_.saturation(oilPhaseIdx),
                        problem.maxOilSaturation(globalSpaceIdx));
        }

        // take the meaning of the switching primary variable into account for the gas
        // and oil phase compositions
        SegmentIndex segIdx_so;
        if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs) {
            const auto& Rs = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
            fluidState_.setRs(Rs);
            saturated[oilPhaseIdx] = false;
        } else {
            if (FluidSystem::enableDissolvedGas()) { // Add So > 0? i.e. if only water set rs = 0)
                OPM_TIMEBLOCK_LOCAL(UpdateSaturatedRs);
                const Evaluation& p = fluidState_.pressure(oilPhaseIdx);
                segIdx_so = oilpvt.saturatedGasDissolutionFactorTable()[pvtRegionIdx].findSegmentIndex(p,/*extrapolate=*/true);
                const Evaluation& RsSat = oilpvt.saturatedGasDissolutionFactorTable()[pvtRegionIdx].eval(p, segIdx_so);
                //const Evaluation& RsSat = enableExtbo ? asImp_().rs() :
                //   FluidSystem::saturatedDissolutionFactor(fluidState_,
                //                                           oilPhaseIdx,
                //                                            pvtRegionIdx,
                //                                            SoMax);
                if(RsMax < RsSat){
                    saturated[oilPhaseIdx] = false;
                }
                fluidState_.setRs(min(RsMax, RsSat));
            }
            else if constexpr (compositionSwitchEnabled){
                fluidState_.setRs(0.0);
                saturated[oilPhaseIdx] = false;
            }
        }
        //saturated[oilPhaseIdx] = true;
        if (saturated[oilPhaseIdx]){
            OPM_TIMEBLOCK_LOCAL(OilSaturatedPvt);
            const Evaluation& p = fluidState_.pressure(oilPhaseIdx);
            // exeption from using saturated tables
            //SegmentIndex segIdx = oilpvt.inverseSaturatedOilBTable()[pvtRegionIdx].findSegmentIndex(p,/*extrapolate=*/true);
            SegmentIndex segIdx = segIdx_so;
            Evaluation b  =oilpvt.inverseSaturatedOilBTable()[pvtRegionIdx].eval(p, segIdx);
            //Evaluation b = oilpvt.saturatedInverseFormationVolumeFactor(pvtRegionIdx, T, p);
            //Evaluation mu = oilpvt.saturatedViscosity(pvtRegionIdx, T, p);
            //Evaluation invBMu = oilpvt.inverseSaturatedOilBMuTable()[pvtRegionIdx].eval(p, /*extrapolate=*/true);
            Evaluation invBMu = oilpvt.inverseSaturatedOilBMuTable()[pvtRegionIdx].eval(p,segIdx);
            Evaluation mu = b/invBMu;
            fluidState_.setInvB(oilPhaseIdx, b);
            viscosity[oilPhaseIdx] =mu;
        }else{
            OPM_TIMEBLOCK_LOCAL(OilUnSaturatedPvt);
            const Evaluation& p = fluidState_.pressure(oilPhaseIdx);
            const Evaluation& Rs = fluidState_.Rs();
            //Evaluation b = oilpvt.inverseFormationVolumeFactor(pvtRegionIdx, T, p,Rs);//??
            unsigned ii,j1,j2;
            Evaluation alpha, beta1, beta2;
            //Evaluation b = oilpvt.inverseOilBTable()[pvtRegionIdx].eval(Rs,p,/*extrapolate*/true);
            oilpvt.inverseOilBTable()[pvtRegionIdx].findPoints(ii,j1,j2,alpha, beta1,beta2,Rs,p,/*extrapolate*/true);
            Evaluation b = oilpvt.inverseOilBTable()[pvtRegionIdx].eval(ii,j1,j2,alpha, beta1,beta2);//,Rs,p,/*extrapolate*/true);
            //Evaluation mu = oilpvt.viscosity(pvtRegionIdx, T, p, Rs);
            //Evaluation invBMu = oilpvt.inverseOilBMuTable()[pvtRegionIdx].eval(Rs, p, /*extrapolate=*/true);
            Evaluation invBMu = oilpvt.inverseOilBMuTable()[pvtRegionIdx].eval(ii,j1,j2,alpha, beta1,beta2);//,Rs,p,/*extrapolate*/true);
            Evaluation mu = b/invBMu;
            fluidState_.setInvB(oilPhaseIdx, b);
            viscosity[oilPhaseIdx] = mu;
        }

        // if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rv) {
        //     const auto& Rv = priVars.makeEvaluation(Indices::compositionSwitchIdx, timeIdx);
        //     fluidState_.setRv(Rv);
        //     saturated[gasPhaseIdx] = false;
        // } else {
        //     if (FluidSystem::enableVaporizedOil() ) { // Add Sg > 0? i.e. if only water set rv = 0)
        //         OPM_TIMEBLOCK_LOCAL(UpdateSaturatedRv);
        //         //NB! should save the indexing for later evalustion
        //         const Evaluation& RvSat = enableExtbo ? asImp_().rv() :
        //             FluidSystem::saturatedDissolutionFactor(fluidState_,
        //                                                     gasPhaseIdx,
        //                                                     pvtRegionIdx,
        //                                                     SoMax);
        //         if(RvMax < RvSat){
        //             saturated[gasPhaseIdx] = false;
        //         }
        //         fluidState_.setRv(min(RvMax, RvSat));

        //     }
        //     else if constexpr (compositionSwitchEnabled){
        //         fluidState_.setRv(0.0);
        //         saturated[gasPhaseIdx] = true;
        //     }
        // }
        //if(saturated[gasPhaseIdx]){
        if(true){
            OPM_TIMEBLOCK_LOCAL(GasSaturatedPvt);
            const Evaluation& p = fluidState_.pressure(gasPhaseIdx);
            // no oil  gas present  and enableVaporized oil
            SegmentIndex segIdx = gaspvt.inverseGasB()[pvtRegionIdx].findSegmentIndex(p,/*extrapolate=*/true);
            Evaluation b  =gaspvt.inverseGasB()[pvtRegionIdx].eval(p, segIdx);
            //Evaluation b = gaspvt.saturatedInverseFormationVolumeFactor(pvtRegionIdx, T, p);
            //Evaluation mu = gaspvt.saturatedViscosity(pvtRegionIdx, T, p);
            //Evaluation invBMu = gaspvt.inverseGasBMu()[pvtRegionIdx].eval(p, /*extrapolate=*/true);
            const Evaluation& invBMu = gaspvt.inverseGasBMu()[pvtRegionIdx].eval(p, segIdx);
            Evaluation mu = b/invBMu;
            fluidState_.setInvB(gasPhaseIdx, b);
            viscosity[gasPhaseIdx] = mu;
        }else{
            OPM_TIMEBLOCK_LOCAL(GasUnSaturatedPvt);
            const Evaluation& p = fluidState_.pressure(gasPhaseIdx);
            const Evaluation& Rv = fluidState_.Rv();
            const Evaluation& Rvw = 0.0;
            Evaluation b = gaspvt.inverseFormationVolumeFactor(pvtRegionIdx, T, p, Rv, Rvw);
            //Evaluation mu = gaspvt.viscosity(pvtRegionIdx, T, p, Rv, Rvw);
            Evaluation invBMu = gaspvt.inverseGasBMu()[pvtRegionIdx].eval(p, /*extrapolate=*/true);
            Evaluation mu = b/invBMu;
            fluidState_.setInvB(gasPhaseIdx, b);
            viscosity[gasPhaseIdx] = mu;
        }




        if (priVars.primaryVarsMeaningWater() == PrimaryVariables::WaterMeaning::Rvw) {
            const auto& Rvw = priVars.makeEvaluation(Indices::waterSwitchIdx, timeIdx);
            fluidState_.setRvw(Rvw);
            saturated[waterPhaseIdx] = false;
        } else {
            //NB! should save the indexing for later evaluation
            if (FluidSystem::enableVaporizedWater()) { // Add Sg > 0? i.e. if only water set rv = 0)
                OPM_TIMEBLOCK_LOCAL(UpdateSaturatedRv);
                const Evaluation& RvwSat = FluidSystem::saturatedVaporizationFactor(fluidState_,
                                                            gasPhaseIdx,
                                                            pvtRegionIdx);
                fluidState_.setRvw(RvwSat);
                saturated[waterPhaseIdx] = true;
            }
        }

        // typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        // paramCache.setRegionIndex(pvtRegionIdx);
        // if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
        //     paramCache.setMaxOilSat(SoMax);
        // }
        // paramCache.updateAll(fluidState_);

        // compute the phase densities and transform the phase permeabilities into mobilities
        // int nmobilities = 1;
        // std::vector<std::array<Evaluation,numPhases>*> mobilities = {&mobility_};
        // if (dirMob_) {
        //     for (int i=0; i<3; i++) {
        //         nmobilities += 1;
        //         mobilities.push_back(&(dirMob_->getArray(i)));
        //     }
        // }

        // for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
        //     if (!FluidSystem::phaseIsActive(phaseIdx))
        //         continue;
        //     computeInverseFormationVolumeFactorAndViscosity(fluidState_, phaseIdx, pvtRegionIdx, SoMax);
        // }

        {
            OPM_TIMEBLOCK_LOCAL(WaterPvt);
            //Evaluation salt= 0.0;
            const Evaluation& p = fluidState_.pressure(waterPhaseIdx);
//            Evaluation b = waterpvt.inverseFormationVolumeFactor(pvtRegionIdx, T, p, salt);
            //, Rsw, saltConcentration);
//            Evaluation mu = waterpvt.viscosity(pvtRegionIdx, T, p, salt);
            //auto& ccWaterPvt = waterpvt.template getRealPvt<Opm::WaterPvtApproach::ConstantCompressibilityWater>();
            Evaluation b;
            Evaluation mu;
            waterpvt.inverseBAndMu(b,mu,pvtRegionIdx, p);
            fluidState_.setInvB(waterPhaseIdx, b);
            viscosity[waterPhaseIdx] = mu;
        }





        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            //if (!FluidSystem::phaseIsActive(phaseIdx))
            //    continue;
            mobility_[phaseIdx] /= viscosity[phaseIdx];
        }


        Valgrind::CheckDefined(mobility_);

        // calculate the phase densities
        Evaluation rho;
        if (FluidSystem::phaseIsActive(waterPhaseIdx)) {
            OPM_TIMEBLOCK_LOCAL(UpdateWDensity);
            rho = fluidState_.invB(waterPhaseIdx);
            rho *= FluidSystem::referenceDensity(waterPhaseIdx, pvtRegionIdx);
            fluidState_.setDensity(waterPhaseIdx, rho);
        }

        if (FluidSystem::phaseIsActive(gasPhaseIdx)) {
            OPM_TIMEBLOCK_LOCAL(UpdateGDensity);
            rho = fluidState_.invB(gasPhaseIdx);
            rho *= FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
            if (FluidSystem::enableVaporizedOil()) {
                rho +=
                    fluidState_.invB(gasPhaseIdx) *
                    fluidState_.Rv() *
                    FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
            }
            if (FluidSystem::enableVaporizedWater()) {
                rho +=
                    fluidState_.invB(gasPhaseIdx) *
                    fluidState_.Rvw() *
                    FluidSystem::referenceDensity(waterPhaseIdx, pvtRegionIdx);
            }
            fluidState_.setDensity(gasPhaseIdx, rho);
        }

        if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
            OPM_TIMEBLOCK_LOCAL(UpdateODensity);
            rho = fluidState_.invB(oilPhaseIdx);
            rho *= FluidSystem::referenceDensity(oilPhaseIdx, pvtRegionIdx);
            if (FluidSystem::enableDissolvedGas()) {
                rho +=
                    fluidState_.invB(oilPhaseIdx) *
                    fluidState_.Rs() *
                    FluidSystem::referenceDensity(gasPhaseIdx, pvtRegionIdx);
            }
            fluidState_.setDensity(oilPhaseIdx, rho);
        }

        // retrieve the porosity from the problem
        referencePorosity_ = problem.porosity(globalSpaceIdx,timeIdx);
        porosity_ = referencePorosity_;
        // the porosity must be modified by the compressibility of the
        // rock...

        Scalar rockCompressibility = problem.rockCompressibility(globalSpaceIdx);
        if (rockCompressibility > 0.0) {
            OPM_TIMEBLOCK_LOCAL(UpdateRockCompressibility);
            Scalar rockRefPressure = problem.rockReferencePressure(globalSpaceIdx);
            Evaluation x;
            if (FluidSystem::phaseIsActive(oilPhaseIdx)) {
                x = rockCompressibility*(fluidState_.pressure(oilPhaseIdx) - rockRefPressure);
            } else if (FluidSystem::phaseIsActive(waterPhaseIdx)){
                x = rockCompressibility*(fluidState_.pressure(waterPhaseIdx) - rockRefPressure);
            } else {
                x = rockCompressibility*(fluidState_.pressure(gasPhaseIdx) - rockRefPressure);
            }
            porosity_ *= 1.0 + x + 0.5*x*x;
        }

        // deal with water induced rock compaction
        porosity_ *= problem.template rockCompPoroMultiplier<Evaluation>(*this, globalSpaceIdx);

        // the MICP processes change the porosity
        // if constexpr (enableMICP){
        //   Evaluation biofilm_ = priVars.makeEvaluation(Indices::biofilmConcentrationIdx, timeIdx, linearizationType);
        //   Evaluation calcite_ = priVars.makeEvaluation(Indices::calciteConcentrationIdx, timeIdx, linearizationType);
        //   porosity_ += - biofilm_ - calcite_;
        // }

        // // deal with salt-precipitation
        // if (enableSaltPrecipitation && priVars.primaryVarsMeaningBrine() == PrimaryVariables::BrineMeaning::Sp) {
        //     Evaluation Sp = priVars.makeEvaluation(Indices::saltConcentrationIdx, timeIdx);
        //     porosity_ *= (1.0 - Sp);
        // }

        rockCompTransMultiplier_ = problem.template rockCompTransMultiplier<Evaluation>(*this, globalSpaceIdx);

        // asImp_().solventPvtUpdate_(elemCtx, dofIdx, timeIdx);
        // asImp_().zPvtUpdate_();
        // asImp_().polymerPropertiesUpdate_(elemCtx, dofIdx, timeIdx);
        // asImp_().updateEnergyQuantities_(elemCtx, dofIdx, timeIdx, paramCache);
        // asImp_().foamPropertiesUpdate_(elemCtx, dofIdx, timeIdx);
        // asImp_().MICPPropertiesUpdate_(elemCtx, dofIdx, timeIdx);
        // asImp_().saltPropertiesUpdate_(elemCtx, dofIdx, timeIdx);

        // update the quantities which are required by the chosen
        // velocity model
//        FluxIntensiveQuantities::update_(elemCtx, dofIdx, timeIdx);

        // update the diffusion specific quantities of the intensive quantities
        // DiffusionIntensiveQuantities::update_(fluidState_, paramCache, *this, timeIdx);
        static_assert(!enableDiffusion);
#ifndef NDEBUG
        // some safety checks in debug mode
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            assert(isfinite(fluidState_.density(phaseIdx)));
            assert(isfinite(fluidState_.saturation(phaseIdx)));
            assert(isfinite(fluidState_.temperature(phaseIdx)));
            assert(isfinite(fluidState_.pressure(phaseIdx)));
            assert(isfinite(fluidState_.invB(phaseIdx)));
        }
        assert(isfinite(fluidState_.Rs()));
        assert(isfinite(fluidState_.Rv()));
#endif
    }


    /*!
     * \copydoc ImmiscibleIntensiveQuantities::fluidState
     */
    const FluidState& fluidState() const
    { return fluidState_; }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::mobility
     */
    const Evaluation& mobility(unsigned phaseIdx) const
    { return mobility_[phaseIdx]; }

    const Evaluation& mobility(unsigned phaseIdx, FaceDir::DirEnum facedir) const
    {
        using Dir = FaceDir::DirEnum;
        if (dirMob_) {
            switch(facedir) {
                case Dir::XPlus:
                    return dirMob_->mobilityX_[phaseIdx];
                case Dir::YPlus:
                    return dirMob_->mobilityY_[phaseIdx];
                case Dir::ZPlus:
                    return dirMob_->mobilityZ_[phaseIdx];
                default:
                    throw std::runtime_error("Unexpected face direction");
            }
        }
        else {
            return mobility_[phaseIdx];
        }

    }

    void computeInverseFormationVolumeFactorAndViscosity(FluidState& fluidState,
                                                         unsigned phaseIdx,
                                                         unsigned pvtRegionIdx,
                                                         const Evaluation& SoMax){
        OPM_TIMEBLOCK_LOCAL(UpdateInverseFormationFactorAndViscosity);
        {
        OPM_TIMEBLOCK_LOCAL(UpdateFormationFactor);
        const auto& b = FluidSystem::inverseFormationVolumeFactor(fluidState, phaseIdx, pvtRegionIdx);
        fluidState_.setInvB(phaseIdx, b);
        }
        {
            OPM_TIMEBLOCK_LOCAL(UpdateViscosity);
        typename FluidSystem::template ParameterCache<Evaluation> paramCache;
        paramCache.setRegionIndex(pvtRegionIdx);
        if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
            paramCache.setMaxOilSat(SoMax);
        }
        paramCache.updateAll(fluidState_);

        const auto& mu = FluidSystem::viscosity(fluidState, paramCache, phaseIdx);
        // for (int i = 0; i<nmobilities; i++) {
        //     if (enableExtbo && phaseIdx == oilPhaseIdx) {
        //         (*mobilities[i])[phaseIdx] /= asImp_().oilViscosity();
        //     }
        //     else if (enableExtbo && phaseIdx == gasPhaseIdx) {
        //         (*mobilities[i])[phaseIdx] /= asImp_().gasViscosity();
        //     }
        //     else {
        //         (*mobilities[i])[phaseIdx] /= mu;
        //     }
        // }
        mobility_[phaseIdx] /=mu;
        }
        //mobility_[phaseIdx] /= 1e-3;
    }

    void computeRelpermAndPC(std::array<Evaluation,numPhases>& mobility,
                             std::array<Evaluation, numPhases>& pC,
                             const Problem& problem,
                             const FluidState& fluidState,
                             unsigned globalSpaceIdx){
        OPM_TIMEBLOCK_LOCAL(UpdateRelperm);
        const auto& materialParams = problem.materialLawParams(globalSpaceIdx);
        MaterialLaw::capillaryPressures(pC, materialParams, fluidState_);
        problem.updateRelperms(mobility, dirMob_, fluidState, globalSpaceIdx);

        //mobility_[waterPhaseIdx] = Sw;
        //mobility_[oilPhaseIdx] = So;
        //mobility_[gasPhaseIdx] = Sg;
    }
    /*!
     * \copydoc ImmiscibleIntensiveQuantities::porosity
     */
    const Evaluation& porosity() const
    { return porosity_; }

    /*!
     * The pressure-dependent transmissibility multiplier due to rock compressibility.
     */
    const Evaluation& rockCompTransMultiplier() const
    { return rockCompTransMultiplier_; }

    /*!
     * \brief Returns the index of the PVT region used to calculate the thermodynamic
     *        quantities.
     *
     * This allows to specify different Pressure-Volume-Temperature (PVT) relations in
     * different parts of the spatial domain. Note that this concept should be seen as a
     * work-around of the fact that the black-oil model does not capture the
     * thermodynamics well enough. (Because there is, err, only a single real world with
     * in which all substances follow the same physical laws and hence the same
     * thermodynamics.) Anyway: Since the ECL file format uses multiple PVT regions, we
     * support it as well in our black-oil model. (Note that, if it is not explicitly
     * specified, the PVT region index is 0.)
     */
    auto pvtRegionIndex() const
        -> decltype(std::declval<FluidState>().pvtRegionIndex())
    { return fluidState_.pvtRegionIndex(); }

    /*!
     * \copydoc ImmiscibleIntensiveQuantities::relativePermeability
     */
    Evaluation relativePermeability(unsigned phaseIdx) const
    {
        // warning: slow
        return fluidState_.viscosity(phaseIdx)*mobility(phaseIdx);
    }

    /*!
     * \brief Returns the porosity of the rock at reference conditions.
     *
     * I.e., the porosity of rock which is not perturbed by pressure and temperature
     * changes.
     */
    Scalar referencePorosity() const
    { return referencePorosity_; }

private:
    friend BlackOilSolventIntensiveQuantities<TypeTag>;
    friend BlackOilExtboIntensiveQuantities<TypeTag>;
    friend BlackOilPolymerIntensiveQuantities<TypeTag>;
    friend BlackOilEnergyIntensiveQuantities<TypeTag>;
    friend BlackOilFoamIntensiveQuantities<TypeTag>;
    friend BlackOilBrineIntensiveQuantities<TypeTag>;
    friend BlackOilMICPIntensiveQuantities<TypeTag>;

    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    FluidState fluidState_;
    Scalar referencePorosity_;
    Evaluation porosity_;
    Evaluation rockCompTransMultiplier_;
    std::array<Evaluation,numPhases> mobility_;

    // Instead of writing a custom copy constructor and a custom assignment operator just to handle
    // the dirMob_ unique ptr member variable when copying BlackOilIntensiveQuantites (see for example
    // updateIntensitiveQuantities_() in fvbaseelementcontext.hh for a copy example) we write the below
    // custom wrapper class CopyablePtr which wraps the unique ptr and makes it copyable.
    //
    // The advantage of this approach is that we avoid having to call all the base class copy constructors and
    // assignment operators explicitly (which is needed when writing the custom copy constructor and assignment
    // operators) which could become a maintenance burden. For example, when adding a new base class (if that should
    // be needed sometime in the future) to BlackOilIntensiveQuantites we could forget to update the copy
    // constructor and assignment operators.
    //
    // We want each copy of the BlackOilIntensiveQuantites to be unique, (TODO: why?) so we have to make a copy
    // of the unique_ptr each time we copy construct or assign to it from another BlackOilIntensiveQuantites.
    // (On the other hand, if a copy could share the ptr with the original, a shared_ptr could be used instead and the
    // wrapper would not be needed)
    DirectionalMobilityPtr dirMob_;
};

} // namespace Opm

#endif

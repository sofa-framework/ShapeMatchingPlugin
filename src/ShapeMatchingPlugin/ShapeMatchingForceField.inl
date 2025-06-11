/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once

#include <sofa/core/behavior/ForceField.inl>
#include <ShapeMatchingPlugin/ShapeMatchingForceField.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/rmath.h>
#include <assert.h>
#include <iostream>

namespace shapematchingplugin
{

template<class DataTypes>
ShapeMatchingForceField<DataTypes>::ShapeMatchingForceField()
: l_rotationFinder(initLink("rotationFinder", "link to the rotation finder"))
, d_stiffness(initData(&d_stiffness, (Real)500, "stiffness", "force stiffness"))
{
}

template<class DataTypes>
void ShapeMatchingForceField<DataTypes>::init()
{
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);

    Inherit::init();

    if (!l_rotationFinder.get())
    {
        sofa::core::sptr< ShapeMatchingRotationFinder<DataTypes>> rotationFinder;
        this->getContext()->get(rotationFinder);
        if (!rotationFinder)
        {
            msg_error() << "did not found any Rotation Finder.";
            return;
        }
        else
        {
            l_rotationFinder.set(rotationFinder);
        }
    }

    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}

template<class DataTypes>
void ShapeMatchingForceField<DataTypes>::addForce(const core::MechanicalParams* /* mparams */ /* PARAMS FIRST */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& /* v */)
{
    if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
        return;

    auto f1 = sofa::helper::getWriteOnlyAccessor(f);
    auto p1 = sofa::helper::getReadAccessor(x);

    type::Mat<3,1,Real> q;
    const auto& vR = l_rotationFinder->getRotations();
    const auto& vcm = l_rotationFinder->getCM();
    const auto& vcm0 = l_rotationFinder->getCM0();
    const auto& vNeighborhood = l_rotationFinder->getNeighborhood();
    const auto& x0 = this->mstate->read(sofa::core::vec_id::read_access::restPosition)->getValue();
    Coord gi, xi;

    f1.resize(p1.size());
    const auto& stiffness = d_stiffness.getValue();
    for (unsigned int i=0; i<vNeighborhood.size(); ++i)
    {
        const Mat3x3& R = vR[i];
        const Coord& Xcm = vcm[i];
        const Coord& Xcm0 = vcm0[i];

        for(const auto& ni : vNeighborhood[i])
        {
            gi = Xcm + R*(x0[ni] - Xcm0);
            xi = p1[ni];
            f1[ni] += (gi - xi) * stiffness;
        }
    }
}

template<class DataTypes>
void ShapeMatchingForceField<DataTypes>::addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& df, const DataVecDeriv& dx)
{
    if (this->d_componentState.getValue() != sofa::core::objectmodel::ComponentState::Valid)
        return;

    auto df1 = sofa::helper::getWriteOnlyAccessor(df);
    auto dp1 = sofa::helper::getReadAccessor(dx);

    const auto& vDR = l_rotationFinder->getDRotations(dx.getValue());
    const auto& vcm0 = l_rotationFinder->getCM0();
    const auto& vNeighborhood = l_rotationFinder->getNeighborhood();
    const auto& x0 = this->mstate->read(sofa::core::vec_id::read_access::restPosition)->getValue();

    df1.resize(dp1.size());
    const Real fact = d_stiffness.getValue()*mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    for (unsigned int i=0; i<vNeighborhood.size(); ++i)
    {
        const Mat3x3& dR = vDR[i];
        const Coord& Xcm0 = vcm0[i];
        Coord dxcm;

        for(const auto& ni : vNeighborhood[i])
        {
            dxcm += dp1[ni];
        }
        
        dxcm /= vNeighborhood[i].size();

        for (const auto& ni : vNeighborhood[i])
        {
            const Coord dgi = dxcm + dR * (x0[ni] - Xcm0); // + R*(x0[ni] - Xcm0);
            const Coord& dxi = dp1[ni];
            df1[ni] += (dgi - dxi)*fact;
        }
    }
}

template<class DataTypes>
SReal ShapeMatchingForceField<DataTypes>::getPotentialEnergy(const core::MechanicalParams* /*mparams*/, const DataVecCoord&  /*x */) const
{
    // NOT IMPLEMENTED
    return 0.0;
}

template<class DataTypes>
void ShapeMatchingForceField<DataTypes>::draw(const core::visual::VisualParams* )
{
}


} // namespace shapematchingplugin

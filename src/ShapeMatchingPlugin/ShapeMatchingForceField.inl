/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#ifndef SOFA_COMPONENT_FORCEFIELD_SHAPEMATCHINGFORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_SHAPEMATCHINGFORCEFIELD_INL

#include <sofa/core/behavior/ForceField.inl>
#include <ShapeMatchingPlugin/ShapeMatchingForceField.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/rmath.h>
#include <assert.h>
#include <iostream>
#include <sofa/simulation/Node.h>

namespace sofa
{

namespace component
{

namespace forcefield
{

template<class DataTypes>
void ShapeMatchingForceField<DataTypes>::addForce(const core::MechanicalParams* /* mparams */ /* PARAMS FIRST */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& /* v */)
{
    sofa::helper::WriteAccessor< core::objectmodel::Data< VecDeriv > > f1 = f;
    sofa::helper::ReadAccessor< core::objectmodel::Data< VecCoord > > p1 = x;

    type::Mat<3,1,Real> q;
    const type::vector<Mat3x3>& vR = rotationFinder->getRotations();
    const VecCoord& vcm = rotationFinder->getCM();
    const VecCoord& vcm0 = rotationFinder->getCM0();
    const typename container::ShapeMatchingRotationFinder<DataTypes>::VecNeighborhood& vNeighborhood = rotationFinder->getNeighborhood();
    const VecCoord& x0 = this->mstate->read(sofa::core::ConstVecCoordId::restPosition())->getValue();
    Coord gi, xi;

    f1.resize(p1.size());

    for (unsigned int i=0; i<vNeighborhood.size(); ++i)
    {
        Mat3x3 R = vR[i];
        Coord Xcm = vcm[i];
        Coord Xcm0 = vcm0[i];
        typename container::ShapeMatchingRotationFinder<DataTypes>::Neighborhood::const_iterator it, itEnd;
        for (it = vNeighborhood[i].begin(), itEnd = vNeighborhood[i].end(); it != itEnd ; ++it)
        {
            unsigned int neighbor_i = *it;
            gi = Xcm + R*(x0[neighbor_i] - Xcm0);
            xi = p1[neighbor_i];
            f1[neighbor_i] += (gi-xi)*stiffness.getValue();
        }
    }
}

template<class DataTypes>
void ShapeMatchingForceField<DataTypes>::addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& df, const DataVecDeriv& dx)
{
    sofa::helper::WriteAccessor< core::objectmodel::Data< VecDeriv > > df1 = df;
    sofa::helper::ReadAccessor< core::objectmodel::Data< VecDeriv > > dp1 = dx;

    const type::vector<Mat3x3>& vDR = rotationFinder->getDRotations(dx.getValue());
    const VecCoord& vcm0 = rotationFinder->getCM0();
    const typename container::ShapeMatchingRotationFinder<DataTypes>::VecNeighborhood& vNeighborhood = rotationFinder->getNeighborhood();
    const VecCoord& x0 = this->mstate->read(sofa::core::ConstVecCoordId::restPosition())->getValue();

    df1.resize(dp1.size());
    const Real fact = stiffness.getValue()*mparams->kFactorIncludingRayleighDamping(this->rayleighStiffness.getValue());
    for (unsigned int i=0; i<vNeighborhood.size(); ++i)
    {
        const Mat3x3& dR = vDR[i];
        Coord Xcm0 = vcm0[i];
        Coord dxcm;

        for (auto it = vNeighborhood[i].cbegin(), itEnd = vNeighborhood[i].cend(); it != itEnd ; ++it)
        {
            unsigned int neighbor_i = *it;
            dxcm += dp1[neighbor_i];
        }
        dxcm /= vNeighborhood[i].size();
        for (auto it = vNeighborhood[i].cbegin(), itEnd = vNeighborhood[i].cend(); it != itEnd ; ++it)
        {
            unsigned int neighbor_i = *it;
            Coord dgi = dxcm + dR * (x0[neighbor_i] - Xcm0); // + R*(x0[neighbor_i] - Xcm0);
            Coord dxi = dp1[neighbor_i];
            df1[neighbor_i] += (dgi - dxi)*fact;
        }
    }
}

template<class DataTypes>
void ShapeMatchingForceField<DataTypes>::draw(const core::visual::VisualParams* )
{
}


} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_SHAPEMATCHINGFORCEFIELD_INL

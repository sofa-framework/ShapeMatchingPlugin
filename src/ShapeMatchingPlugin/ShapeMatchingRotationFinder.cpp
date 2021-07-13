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
/*
 * RotationFinder.cpp
 *
 *  Created on: 14 avr. 2009
 *      Author: froy
 */
#define SOFA_COMPONENT_CONTAINER_SHAPEMATCHINGROTATIONFINDER_CPP
#include <sofa/defaulttype/Mat.h>
#include <sofa/component/initMiscForcefieldDev.h>
#include <sofa/component/container/ShapeMatchingRotationFinder.inl>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace container
{

using namespace sofa::defaulttype;

SOFA_DECL_CLASS(ShapeMatchingRotationFinder)

// Register in the Factory
int ShapeMatchingRotationFinderClass = core::RegisterObject("ShapeMatchingRotationFinder")
#ifndef SOFA_FLOAT
.add< ShapeMatchingRotationFinder< Vec3dTypes > >()
//.add< ShapeMatchingRotationFinder< Vec2dTypes > >()
#endif
#ifndef SOFA_DOUBLE
.add< ShapeMatchingRotationFinder< Vec3fTypes > >()
//.add< ShapeMatchingRotationFinder< Vec2fTypes > >()
#endif
.addAlias("RotationFinder")
;
#ifndef SOFA_FLOAT
template class SOFA_MISC_FORCEFIELD_DEV_API ShapeMatchingRotationFinder< defaulttype::Vec3dTypes >;
//template class SOFA_MISC_FORCEFIELD_DEV_API ShapeMatchingRotationFinder< defaulttype::Vec2dTypes >;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_MISC_FORCEFIELD_DEV_API ShapeMatchingRotationFinder< defaulttype::Vec3fTypes >;
//template class SOFA_MISC_FORCEFIELD_DEV_API ShapeMatchingRotationFinder< defaulttype::Vec2fTypes >;
#endif

} // namespace container

} // namespace component

} // namespace sofa
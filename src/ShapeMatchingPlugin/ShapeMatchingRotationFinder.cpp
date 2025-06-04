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
/*
 * RotationFinder.cpp
 *
 *  Created on: 14 avr. 2009
 *      Author: froy
 */
#define SOFA_COMPONENT_CONTAINER_SHAPEMATCHINGROTATIONFINDER_CPP
#include <ShapeMatchingPlugin/config.h>

#include <sofa/type/Mat.h>
#include <ShapeMatchingPlugin/ShapeMatchingRotationFinder.inl>
#include <sofa/core/ObjectFactory.h>

namespace shapematchingplugin
{

using namespace sofa::defaulttype;

// Register in the Factory
void registerShapeMatchingRotationFinder(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("ShapeMatchingRotationFinder")
    .add< ShapeMatchingRotationFinder< Vec3Types > >()
    .addAlias("RotationFinder"));
}

template class SOFA_SHAPEMATCHINGPLUGIN_API ShapeMatchingRotationFinder< defaulttype::Vec3Types >;

} // namespace shapematchingplugin

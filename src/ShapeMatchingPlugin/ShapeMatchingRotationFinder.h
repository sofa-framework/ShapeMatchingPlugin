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
 * RotationFinder.h
 *
 *  Created on: 14 avr. 2009
 *      Author: froy
 */

#ifndef SOFA_COMPONENT_CONTAINER_SHAPEMATCHINGROTATIONFINDER_H
#define SOFA_COMPONENT_CONTAINER_SHAPEMATCHINGROTATIONFINDER_H

#include <ShapeMatchingPlugin/config.h>
#include <sofa/core/behavior/RotationFinder.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <sofa/defaulttype/VecTypes.h>

/*
 *	This class find Rotation Matrix from two position states (rest and current state)
 *	got from an associated Mechanical State.
 *	Formula & method taken from :
 *	"Meshless Deformations Based on Shape Matching" (Muller, Heidelberger, Teschner and Gross)
 *
 */

namespace sofa
{

namespace core {
	namespace objectmodel {
		class Event;
	} // namespace objectmodel
} // namespace core

namespace component
{

namespace container
{

template <class DataTypes>
class ShapeMatchingRotationFinder : public sofa::core::behavior::RotationFinder<DataTypes>
{
public:
	SOFA_CLASS(ShapeMatchingRotationFinder, SOFA_TEMPLATE(sofa::core::behavior::RotationFinder, DataTypes));

	typedef typename DataTypes::VecCoord VecCoord;
	typedef typename DataTypes::VecDeriv VecDeriv;
	typedef typename DataTypes::Coord Coord;
	typedef typename DataTypes::Deriv Deriv;
	typedef typename Coord::value_type Real;
	typedef type::Mat<3,3,Real> Mat3x3;
	typedef type::fixed_array<Mat3x3,9> DMat3x3;

	typedef core::topology::BaseMeshTopology::PointID Point;

	typedef std::set<Point> Neighborhood;
	typedef type::vector<Neighborhood> VecNeighborhood;

private:
    core::behavior::MechanicalState<DataTypes>* mechanicalState;
    core::topology::BaseMeshTopology* topo;

	//rest data
	Coord x0_cm;
	unsigned int oldRestPositionSize;
	VecNeighborhood pointNeighborhood, lastPointNeighborhood;
	VecCoord Xcm, Xcm0;
	type::vector<Mat3x3> rotations, vecA;
	type::vector<DMat3x3> dRotations;
	type::vector<Mat3x3> dRotations_dx;

    template<class T>
    T min(const T a, const T b) const
	{
		return a < b ? a : b;
	}

	void ComputeNeighborhoodFromNeighborhood();

public:
	enum Axis { X, Y, Z };

	Data<int> axisToFlip;
	Data<bool> showRotations;
	Data<int> neighborhoodLevel;
	Data<int> numOfClusters;
	Data<unsigned int> maxIter;
	Data<Real> epsilon;
	Data<Real> radius;

protected:
	ShapeMatchingRotationFinder();
	virtual ~ShapeMatchingRotationFinder();

public:
	void init();

	const type::vector<Mat3x3>& getRotations();
	
	void getRotations(defaulttype::BaseMatrix * m,int offset = 0) ;

	const type::vector<DMat3x3>& getDRotations();

	const type::vector<Mat3x3>& getDRotations(const VecDeriv& dx);

	void computeQT();

	const VecCoord& getCM() { return Xcm; }

	const VecCoord& getCM0() { return Xcm0; }

	void computeNeighborhood();

	const VecNeighborhood& getNeighborhood();

	void handleEvent(sofa::core::objectmodel::Event*);

	void draw(const core::visual::VisualParams* vparams);

	void flipAxis(Mat3x3 & rotation);
        
public:
	/// Pre-construction check method called by ObjectFactory.
	/// Check that DataTypes matches the MechanicalState.
	template<class T>
	static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
	{
		if (dynamic_cast< core::behavior::MechanicalState< DataTypes >* >(context->getMechanicalState()) == NULL)
			return false;

		return core::objectmodel::BaseObject::canCreate(obj, context, arg);
	}

	virtual std::string getTemplateName() const
	{
		return templateName(this);
	}

	static std::string templateName(const ShapeMatchingRotationFinder<DataTypes>* = NULL)
	{
		return DataTypes::Name();
	}
};

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_CONTAINER_SHAPEMATCHINGROTATIONFINDER_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_SHAPEMATCHINGPLUGIN_API ShapeMatchingRotationFinder< defaulttype::Vec3dTypes >;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_SHAPEMATCHINGPLUGIN_API ShapeMatchingRotationFinder< defaulttype::Vec3fTypes >;
#endif
#endif

} // namespace container

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_CONTAINER_SHAPEMATCHINGROTATIONFINDER_H */

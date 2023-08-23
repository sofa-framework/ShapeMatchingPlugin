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

#include <ShapeMatchingPlugin/config.h>
#include <sofa/core/behavior/RotationFinder.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <sofa/defaulttype/VecTypes.h>
#include <unordered_set>

/*
 *	This class find Rotation Matrix from two position states (rest and current state)
 *	got from an associated Mechanical State.
 *	Formula & method taken from :
 *	"Meshless Deformations Based on Shape Matching" (Muller, Heidelberger, Teschner and Gross)
 *
 */
namespace sofa::component::container
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

	typedef std::unordered_set<Point> Neighborhood;
	typedef type::vector<Neighborhood> VecNeighborhood;

private:
	SingleLink<ShapeMatchingRotationFinder<DataTypes>, core::behavior::MechanicalState<DataTypes>, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_mechanicalState;
	SingleLink<ShapeMatchingRotationFinder<DataTypes>, sofa::core::topology::BaseMeshTopology, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_topology;
    
	//rest data
	std::size_t m_oldRestPositionSize;
	VecNeighborhood m_pointNeighborhood, m_lastPointNeighborhood;
	VecCoord m_Xcm, m_Xcm0;
	type::vector<Mat3x3> m_rotations, m_vecA;
	type::vector<DMat3x3> m_dRotations;
	type::vector<Mat3x3> m_dRotations_dx;

	void ComputeNeighborhoodFromNeighborhood();

public:
	enum class Axis { X, Y, Z };

	Data<int> d_axisToFlip;
	Data<bool> d_showRotations;
	Data<int> d_neighborhoodLevel;
	Data<unsigned> d_numOfClusters;
	Data<unsigned> d_maxIter;
	Data<Real> d_epsilon;
	Data<Real> d_radius;

protected:
	ShapeMatchingRotationFinder();
	~ShapeMatchingRotationFinder() override;

public:
	void init() override;
	void draw(const core::visual::VisualParams* vparams) override;
	void getRotations(linearalgebra::BaseMatrix* m, int offset = 0) override;

	const type::vector<Mat3x3>& getRotations() override;
	const type::vector<DMat3x3>& getDRotations();
	const type::vector<Mat3x3>& getDRotations(const VecDeriv& dx);

	void computeQT();
	const VecCoord& getCM() { return m_Xcm; }
	const VecCoord& getCM0() { return m_Xcm0; }
	void computeNeighborhood();
	const VecNeighborhood& getNeighborhood();
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
};

#if !defined(SOFA_COMPONENT_CONTAINER_SHAPEMATCHINGROTATIONFINDER_CPP)
extern template class SOFA_SHAPEMATCHINGPLUGIN_API ShapeMatchingRotationFinder< defaulttype::Vec3Types >;
#endif

} // namespace sofa::component::container

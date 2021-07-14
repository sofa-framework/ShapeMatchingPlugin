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
#pragma once

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/core/objectmodel/Event.h>
#include <ShapeMatchingPlugin/ShapeMatchingRotationFinder.h>

namespace sofa::component::forcefield
{

/// Meshless deformations based on shape matching
/// https://www.researchgate.net/publication/220184721_Meshless_deformations_based_on_shape_matching
/// 
template<class DataTypes>
class ShapeMatchingForceField : public core::behavior::ForceField<DataTypes>
{
public:
  SOFA_CLASS(SOFA_TEMPLATE(ShapeMatchingForceField, DataTypes), SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

	typedef core::behavior::ForceField<DataTypes> Inherit;
	typedef typename DataTypes::VecCoord VecCoord;
	typedef typename DataTypes::VecDeriv VecDeriv;
	typedef typename DataTypes::Coord Coord;
	typedef typename DataTypes::Deriv Deriv;
	typedef typename Coord::value_type Real;
	typedef type::Mat<3,3,Real> Mat3x3;

    typedef core::objectmodel::Data<VecCoord> DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;
	
protected:
    ShapeMatchingForceField();

public:
    Data<Real> d_stiffness;
    container::ShapeMatchingRotationFinder<DataTypes>* m_rotationFinder;

    virtual void addForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v) override;
    virtual void addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& df, const DataVecDeriv& dx) override;
    virtual SReal getPotentialEnergy(const core::MechanicalParams* /*mparams*/, const DataVecCoord&  /*x */) const override;

    void draw(const core::visual::VisualParams* vparams) override;

    /// Pre-construction check method called by ObjectFactory.
    template<class T>
    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
		using namespace container;

        if (dynamic_cast< ShapeMatchingRotationFinder<DataTypes>* >(context->get< ShapeMatchingRotationFinder<DataTypes> >()) == NULL)
            return false;

        return core::objectmodel::BaseObject::canCreate(obj, context, arg);
    }

    /// Construction method called by ObjectFactory.
    template<class T>
	static typename T::SPtr create(T* tObj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg )
    {
		using namespace container;

        typename T::SPtr obj = core::objectmodel::BaseObject::create(tObj, context, arg);

        if (context)
        {
            obj->m_rotationFinder = dynamic_cast< ShapeMatchingRotationFinder<DataTypes>* >(context->get< ShapeMatchingRotationFinder<DataTypes> >());
        }


		return obj;
    }

    virtual std::string getTemplateName() const
      {
        return templateName(this);
      }

    static std::string templateName(const ShapeMatchingForceField<DataTypes>* = NULL)
    {
      return DataTypes::Name();
    }

};

#if !defined(SOFA_COMPONENT_FORCEFIELD_SHAPEMATCHINGFORCEFIELD_CPP)
extern template class SOFA_SHAPEMATCHINGPLUGIN_API ShapeMatchingForceField<defaulttype::Vec3dTypes>;
#endif // !defined(SOFA_COMPONENT_FORCEFIELD_SHAPEMATCHINGFORCEFIELD_CPP)

} // namespace sofa::component::forcefield

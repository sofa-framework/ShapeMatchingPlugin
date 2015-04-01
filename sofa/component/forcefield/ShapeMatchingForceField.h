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
#ifndef SOFA_COMPONENT_FORCEFIELD_SHAPEMATCHINGFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_SHAPEMATCHINGFORCEFIELD_H

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/component/container/ShapeMatchingRotationFinder.h>

#include <sofa/component/component.h>

namespace sofa
{

namespace component
{

namespace forcefield
{

/// This class can be overridden if needed for additionnal storage within template specializations.
template<class DataTypes>
class ShapeMatchingForceFieldInternalData
{
public:
};

/// Meshless deformations based on shape matching

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
	typedef defaulttype::Mat<3,3,Real> Mat3x3;

    typedef core::objectmodel::Data<VecCoord> DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;
	
protected:
    class Contact
    {
    public:
        int index;
        Coord normal;
        Real fact;
	Contact( int index=0, Coord normal=Coord(),Real fact=Real(0))
	  : index(index),normal(normal),fact(fact)
	  {
	  }

        inline friend std::istream& operator >> ( std::istream& in, Contact& c ){
	  in>>c.index>>c.normal>>c.fact;
            return in;
        }

        inline friend std::ostream& operator << ( std::ostream& out, const Contact& c ){
	  out << c.index << " " << c.normal << " " << c.fact ;
	  return out;
        }

    };

    Data<sofa::helper::vector<Contact> > contacts;

    core::behavior::MechanicalState<DataTypes> * centerDOF;

    ShapeMatchingForceFieldInternalData<DataTypes> data;

	ShapeMatchingForceField()
      : contacts(initData(&contacts,"contacts", "Contacts"))
      , centerDOF(NULL)
      , sphereCenter(initData(&sphereCenter, "center", "sphere center"))
      , sphereRadius(initData(&sphereRadius, (Real)1, "radius", "sphere radius"))
      , stiffness(initData(&stiffness, (Real)500, "stiffness", "force stiffness"))
      , damping(initData(&damping, (Real)5, "damping", "force damping"))
      , color(initData(&color, defaulttype::Vec3f(0.0f,0.0f,1.0f), "color", "sphere color"))
      , bDraw(initData(&bDraw, true, "draw", "enable/disable drawing of the sphere"))
      , centerState(initData(&centerState, "centerState", "path to the MechanicalState controlling the center point"))
      , filter(initData(&filter, (Real)0, "filter", "filter"))
    {
    }

public:

    Data<Coord> sphereCenter;
    Data<Real> sphereRadius;
    Data<Real> stiffness;
    Data<Real> damping;
    Data<defaulttype::Vec3f> color;
    Data<bool> bDraw;
    Data<std::string> centerState;
    Data < Real > filter;
    container::ShapeMatchingRotationFinder<DataTypes>* rotationFinder;

    void setStiffness(Real stiff)
    {
        stiffness.setValue( stiff );
    }

    void setDamping(Real damp)
    {
        damping.setValue( damp );
    }

    virtual void addForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v);
    virtual void addDForce(const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& df, const DataVecDeriv& dx);

    virtual SReal getPotentialEnergy(const core::MechanicalParams* /*mparams*/, const DataVecCoord&  /*x */) const
    {
        // NOT IMPLEMENTED
        return 0.0;
    }

    void draw(const core::visual::VisualParams* vparams);

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
            obj->rotationFinder = dynamic_cast< ShapeMatchingRotationFinder<DataTypes>* >(context->get< ShapeMatchingRotationFinder<DataTypes> >());
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

using sofa::defaulttype::Vec3dTypes;
using sofa::defaulttype::Vec3fTypes;

#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_SHAPEMATCHINGFORCEFIELD_CPP)

#ifndef SOFA_FLOAT
extern template class SOFA_MISC_FORCEFIELD_DEV_API ShapeMatchingForceField<Vec3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_MISC_FORCEFIELD_DEV_API ShapeMatchingForceField<Vec3fTypes>;
#endif

#endif // defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_SHAPEMATCHINGFORCEFIELD_CPP)

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_SHAPEMATCHINGFORCEFIELD_H

/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
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
//
// C++ Interface: EvolutionParameterController
//
// Description:
//
//
// Author: Alex
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef SOFA_COMPONENT_CONTROLLER_AWESOMECONTROLLER_H
#define SOFA_COMPONENT_CONTROLLER_AWESOMECONTROLLER_H

#include <sofa/component/controller/Controller.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/component/component.h>

namespace sofa
{

namespace component
{

namespace controller
{
	
using namespace sofa::defaulttype;

/**
 * @brief EvolutionParameterController Class
 * 
 */
template<class DataTypes>
class EvolutionParameterController : public Controller
{
public:
  SOFA_CLASS(SOFA_TEMPLATE(EvolutionParameterController,DataTypes),Controller);
	typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;

    Data<int> f_startDTAppl;
    Data<int> f_numDTAppl;
    Data<Real> f_controlFactor, f_compControlFactor;
    int simulationStep;

	/**
	 * @brief Default Constructor.
	 */
    EvolutionParameterController();

	/**
	 * @brief Default Destructor.
	 */
    virtual ~EvolutionParameterController(){}

	/**
	 * @brief SceneGraph callback initialization method.
	 */
	void init();

	/**
	 * @name Controller Interface
	 */
	//@{

	/**
	 * @brief Begin Animation event callback.
	 */
	void onBeginAnimationStep(const double dt);

	//@}

	/**
	 * @name Accessors
	 */
	//@{


	//@}

	/**
	 * @brief Apply the controller modifications to the controlled MechanicalState.
	 */
    void applyController(void);

    virtual std::string getTemplateName() const
    {
      return templateName(this);
    }

    static std::string templateName(const EvolutionParameterController<DataTypes>* = NULL)
    {
      return DataTypes::Name();
    }

protected:

};

#if defined(WIN32) && !defined(SOFA_COMPONENT_CONTROLLER_AWESOMECONTROLLER_CPP)
#pragma warning(disable : 4231)
#ifndef SOFA_FLOAT
extern template class SOFA_COMPONENT_CONTROLLER_API EvolutionParameterController<defaulttype::Vec3dTypes>;
extern template class SOFA_COMPONENT_CONTROLLER_API EvolutionParameterController<defaulttype::Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_COMPONENT_CONTROLLER_API EvolutionParameterController<defaulttype::Vec3fTypes>;
extern template class SOFA_COMPONENT_CONTROLLER_API EvolutionParameterController<defaulttype::Rigid3fTypes>;
#endif
#endif

} // namespace controller

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONTROLLER_AWESOMECONTROLLER_H

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
//
// C++ Interface: MappedBeamToTetraForceField
//
// Description: 
//
//
// Author: The SOFA team </www.sofa-framework.org>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_COMPONENT_FORCEFIELD_MAPPEDBEAMTOTETRAFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_MAPPEDBEAMTOTETRAFORCEFIELD_H

#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/component/forcefield/PlaneForceField.h>
#include <vector>

#include <sofa/component/forcefield/BeamFEMForceField.h>
#include <sofa/component/mapping/BarycentricMapping.h>

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/component/linearsolver/CompressedRowSparseMatrix.h>
#include <sofa/component/linearsolver/DefaultMultiMatrixAccessor.h>

namespace sofa
{

namespace component
{

namespace forcefield
{

using namespace sofa::core;
using namespace helper::system::thread;

/// A box of 6 PlaneForceField that can rotate
template<class DataTypes>
class MappedBeamToTetraForceField : public core::behavior::ForceField<DataTypes>
{
public:
	SOFA_CLASS(SOFA_TEMPLATE(MappedBeamToTetraForceField, DataTypes), SOFA_TEMPLATE(core::behavior::ForceField, DataTypes));

	typedef core::behavior::ForceField<DataTypes> Inherit;
	typedef typename DataTypes::VecCoord VecCoord;
	typedef typename DataTypes::VecDeriv VecDeriv;
	typedef typename DataTypes::Coord Coord;
	typedef typename DataTypes::Deriv Deriv;
	typedef Data<VecCoord>  DataVecCoord;
	typedef Data<VecDeriv>  DataVecDeriv;
	typedef typename Coord::value_type Real;
	typedef typename defaulttype::Rigid3dTypes Rigid;
        typedef sofa::component::mapping::BarycentricMapping< DataTypes, Rigid3dTypes >   BarycentricMapping3_VtoR;

	typedef defaulttype::Mat<6,3,Real> MappingBlockType;
	typedef defaulttype::Mat<6,6,Real> BeamBlockType;
	typedef defaulttype::Mat<3,6,Real> TempBlockType;
	typedef defaulttype::Mat<3,3,Real> TetraBlockType;

	typedef typename sofa::component::linearsolver::CompressedRowSparseMatrix<MappingBlockType>::ColBlockConstIterator MappingColBlockConstIterator;
	typedef typename sofa::component::linearsolver::CompressedRowSparseMatrix<BeamBlockType>::ColBlockConstIterator BeamColBlockConstIterator;
	typedef typename sofa::component::linearsolver::CompressedRowSparseMatrix<TempBlockType>::ColBlockConstIterator TempColBlockConstIterator;

	typedef typename sofa::component::linearsolver::CompressedRowSparseMatrix<MappingBlockType>::RowBlockConstIterator MappingRowBlockConstIterator;
	typedef typename sofa::component::linearsolver::CompressedRowSparseMatrix<BeamBlockType>::RowBlockConstIterator BeamRowBlockConstIterator;
	typedef typename sofa::component::linearsolver::CompressedRowSparseMatrix<TempBlockType>::RowBlockConstIterator TempRowBlockConstIterator;

	typedef typename sofa::component::linearsolver::CompressedRowSparseMatrix<MappingBlockType>::BlockConstAccessor MappingBlockConstAccessor;
	typedef typename sofa::component::linearsolver::CompressedRowSparseMatrix<BeamBlockType>::BlockConstAccessor BeamBlockConstAccessor;
	typedef typename sofa::component::linearsolver::CompressedRowSparseMatrix<TempBlockType>::BlockConstAccessor TempBlockConstAccessor;

	typedef typename sofa::component::linearsolver::CompressedRowSparseMatrix<MappingBlockType>::BlockAccessor MappingBlockAccessor;
	typedef typename sofa::component::linearsolver::CompressedRowSparseMatrix<BeamBlockType>::BlockAccessor BeamBlockAccessor;
	typedef typename sofa::component::linearsolver::CompressedRowSparseMatrix<TempBlockType>::BlockAccessor TempBlockAccessor;

        Data<std::string> _beamMO;
        Data<std::string> _beamBaryMapping;
        Data<std::string> _beamFEM;        

protected:
	BeamFEMForceField<Rigid> *mappedBeamForceField;
	BarycentricMapping3_VtoR *mappingBeamTetra;
	core::behavior::MechanicalState<DataTypes>* parenchymaMO;
	core::behavior::MechanicalState<Rigid>* beamMO;
        sofa::component::linearsolver::DefaultMultiMatrixAccessor* beamMatrixAccessor;
        sofa::component::linearsolver::CompressedRowSparseMatrix<BeamBlockType>* beamMatrix ;
        unsigned int nbBeamDOFs, nbTetDOFs;
        sofa::helper::system::thread::CTime *timer;

	MappedBeamToTetraForceField(core::behavior::MechanicalState<DataTypes>* /*object*/=NULL, const std::string& /*name*/="")               
        : _beamMO( initData(&_beamMO, "beamMechObject", "beam mechanical object mapped to the tetra") )
        , _beamBaryMapping( initData(&_beamBaryMapping, "beamBarycentricMapping", "barycentric mapping between the beam and tetra") )
        , _beamFEM( initData(&_beamFEM, "beamFEM", "beam FEM force field") )
	{		

	}        

	~MappedBeamToTetraForceField()
	{
	}

public:
	virtual void bwdInit();


	virtual void addForce (const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v);

	virtual void addDForce (const core::MechanicalParams* mparams /* PARAMS FIRST */, DataVecDeriv& df, const DataVecDeriv& dx);

	virtual void addKToMatrix(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, const sofa::core::behavior::MultiMatrixAccessor* matrix );

	virtual double getPotentialEnergy(const MechanicalParams* /*mparams*/ /* PARAMS FIRST */, const DataVecCoord&   x) const;


//	virtual void addForce (VecDeriv& f, const VecCoord& x, const VecDeriv& v);
//
//	virtual void addDForce (VecDeriv& df, const VecDeriv& dx, double kFactor, double bFactor);
//
//	virtual double getPotentialEnergy(const VecCoord& x) const;
//
//	void addKToMatrix(sofa::defaulttype::BaseMatrix *mat, SReal k, unsigned int &offset);

	void draw(const core::visual::VisualParams* vparams);

};

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif

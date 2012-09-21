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
#ifndef SOFA_COMPONENT_FORCEFIELD_MAPPEDBEAMTOTETRAFORCEFIELD_INL
#define SOFA_COMPONENT_FORCEFIELD_MAPPEDBEAMTOTETRAFORCEFIELD_INL

#include <sofa/core/behavior/ForceField.inl>
#include <sofa/helper/system/config.h>
#include <sofa/helper/system/gl.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sofa/component/forcefield/MappedBeamToTetraForceField.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/system/thread/CTime.h>

namespace sofa
{

namespace component
{

namespace forcefield
{
template<class DataTypes>
void MappedBeamToTetraForceField<DataTypes>::bwdInit()
{
    this->core::behavior::ForceField<DataTypes>::bwdInit();
    sout << "Initializing Mapped BeamToTetra Force Field" << sendl;

    this->getContext()->get(mappedBeamForceField, _beamFEM.getValue());
    if (mappedBeamForceField == NULL) {
        serr << "No BeamForceField  to mapped found " << sendl;
    } else {
        sout << "Beam force field found, name = " << mappedBeamForceField->getName() << sendl;
    }

    this->getContext()->get(mappingBeamTetra, _beamBaryMapping.getValue());
    if (mappingBeamTetra == NULL) {
        serr << "No Mapping found" << sendl;
    } else {
        sout <<  "Barycentric mapping found: " << mappingBeamTetra->getName() << sendl;
    }

    this->getContext()->get(beamMO, _beamMO.getValue());
    if (!beamMO)
        serr << "Beam mechanical object not found! " << sendl;
    else
        sout << "Beam mechanical object found: " << beamMO->getName() << sendl;

    nbBeamDOFs = beamMO->getSize();
    nbTetDOFs =  this->getMState()->getSize();

    beamMatrix = new sofa::component::linearsolver::CompressedRowSparseMatrix<BeamBlockType>;
    std::cout << "Size of beam matrix: " << nbBeamDOFs << std::endl;
    beamMatrix->resizeBloc(nbBeamDOFs,nbBeamDOFs);

    beamMatrixAccessor= new sofa::component::linearsolver::DefaultMultiMatrixAccessor;
    beamMatrixAccessor->addMechanicalState(beamMO);
    beamMatrixAccessor->setGlobalMatrix(beamMatrix);
    beamMatrixAccessor->setupMatrices();

    timer = new CTime();

}

template<class DataTypes>
void MappedBeamToTetraForceField<DataTypes>::addForce(const core::MechanicalParams*  /* mparams*/, DataVecDeriv& , const DataVecCoord& , const DataVecDeriv& )
{    
    //std::cout << "AddForce in MappedFF" << std::endl;
}

template<class DataTypes>
void MappedBeamToTetraForceField<DataTypes>::addDForce(const core::MechanicalParams*  /* mparams */, DataVecDeriv& , const DataVecDeriv& )
{    
    //double k = mparams->kFactor();
    //std::cout << "AddDForce in MappedFF" << std::endl;
}

////void MappedBeamToTetraForceField<DataTypes>::addKToMatrix(sofa::defaulttype::BaseMatrix *mat, SReal k, unsigned int &offset)
template<class DataTypes>
void MappedBeamToTetraForceField<DataTypes>::addKToMatrix(const sofa::core::MechanicalParams* mparams /* PARAMS FIRST */, const sofa::core::behavior::MultiMatrixAccessor* matrix )
{
    //double startTime, stopTime;
    //startTime = (double)timer->getTime();
    sofa::core::behavior::MultiMatrixAccessor::MatrixRef r = matrix->getMatrix(this->mstate);
    sofa::defaulttype::BaseMatrix* mat = r.matrix;

    beamMatrix->clear();
    mappedBeamForceField->addKToMatrix(mparams,beamMatrixAccessor);

    const sofa::component::linearsolver::CompressedRowSparseMatrix<MappingBlockType>* mappingMatrix;
    mappingMatrix = dynamic_cast<const sofa::component::linearsolver::CompressedRowSparseMatrix<MappingBlockType> *>(mappingBeamTetra->getJ());

    MappingBlockType mappingBuffer;
    BeamBlockType beamBuffer;

    sofa::component::linearsolver::CompressedRowSparseMatrix<TempBlockType> tempMatrix;
    tempMatrix.resizeBloc(nbTetDOFs, nbBeamDOFs);

    //stopTime = (double)timer->getTime();
    //std::cout << "Length1 = " << stopTime - startTime << std::endl;
    //startTime = (double)timer->getTime();

    /*
     * K11 = mat
     * K22 = beamMatrix
     *   J = mappingMatrix
     * */
    //operation K11 += Jt * K22 * J evaluated in two steps

    //op1 : tempMatrix = Jt * K22
    for (int mappingRowIndex = 0; mappingRowIndex < mappingMatrix->nBlocRow; mappingRowIndex++) {   //through X, must take each row  (but test can be added)
        for (MappingColBlockConstIterator mappingColIter = mappingMatrix->bRowBegin(mappingRowIndex); mappingColIter < mappingMatrix->bRowEnd(mappingRowIndex); mappingColIter++) {  //take non zero blocks in row, determines the row in K)
            MappingBlockConstAccessor mappingBlock = mappingColIter.bloc();
            const MappingBlockType& mappingBlockData = *(const MappingBlockType*)mappingBlock.elements(mappingBuffer.ptr());
            int mappingColIndex = mappingBlock.getCol();
            for (BeamColBlockConstIterator beamColIter = beamMatrix->bRowBegin(mappingRowIndex); beamColIter < beamMatrix->bRowEnd(mappingRowIndex); beamColIter++) {
                BeamBlockConstAccessor beamBlock = beamColIter.bloc();
                const BeamBlockType& beamBlockData = *(const BeamBlockType*)beamBlock.elements(beamBuffer.ptr());
                int beamColIndex = beamBlock.getCol();
                TempBlockType tempBlockData(0.0);
                //multiply the block, could be done more efficiently
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 6; j++)
                        for (int k = 0; k < 6; k++)
                            tempBlockData(i,j) += mappingBlockData(k,i)*beamBlockData(k,j);
                tempMatrix.blocAdd(mappingColIndex, beamColIndex, tempBlockData.ptr());
            }
        }
    }

    TempBlockType tempBuffer;
    //op1 : K11 += tempMatrix * J
    for (int tempRowIndex = 0; tempRowIndex < tempMatrix.nBlocRow; tempRowIndex++) {
        for (TempColBlockConstIterator tempColIter = tempMatrix.bRowBegin(tempRowIndex); tempColIter < tempMatrix.bRowEnd(tempRowIndex); tempColIter++) {
            TempBlockConstAccessor tempBlock = tempColIter.bloc();
            const TempBlockType& tempBlockData = *(const TempBlockType*) tempBlock.elements(tempBuffer.ptr());
            int tempColIndex = tempBlock.getCol();
            for (MappingColBlockConstIterator mappingColIter = mappingMatrix->bRowBegin(tempColIndex); mappingColIter < mappingMatrix->bRowEnd(tempColIndex); mappingColIter++) {
                MappingBlockConstAccessor mappingBlock = mappingColIter.bloc();
                const MappingBlockType &mappingBlockData = *(const MappingBlockType*) mappingBlock.elements(mappingBuffer.ptr());
                int mappingColIndex = mappingBlock.getCol();
                TetraBlockType tetraBlockData(0.0);
                //multiply the block, could be done more efficiently
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++) {
                        for (int k = 0; k < 6; k++)
                            tetraBlockData(i,j) += tempBlockData(i,k) * mappingBlockData(k,j);
                        //because input matrix is  BaseMatrix, we do it like this for now
                        mat->add(3*tempRowIndex+i, 3*mappingColIndex+j, tetraBlockData(i,j));
                    }
                //mat->blocAdd(tempRowIndex,mappingColIndex,tetraBlockData.ptr());   //if mat is block matrix
            }
        }
    }
    //stopTime = (double)timer->getTime();
    //std::cout << "Length = " << stopTime - startTime << std::endl;
}


template <class DataTypes>
double MappedBeamToTetraForceField<DataTypes>::getPotentialEnergy(const MechanicalParams* /* PARAMS FIRST */, const DataVecCoord& ) const
{
    return(0.0);
}

template<class DataTypes>
void MappedBeamToTetraForceField<DataTypes>::draw(const core::visual::VisualParams* )
{
}

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif

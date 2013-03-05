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

    this->getContext()->get(mappedForceField, _mappedFEM.getValue());
    if (mappedForceField == NULL) {
        serr << "No mapped ForceField found " << sendl;
    } else {
        sout << "Found force field, name = " << mappedForceField->getName() << sendl;
    }

    this->getContext()->get(mapping, _mapping.getValue());
    if (mapping == NULL) {
        serr << "No Mapping found" << sendl;
    } else {
        sout <<  "Mapping found: " << mapping->getName() << sendl;
    }

    this->getContext()->get(mappedMO, _mappedMO.getValue());
    if (!mappedMO)
        serr << "Mapped Mechanical object not found! " << sendl;
    else
        sout << "Mapped Mechanical object found: " << mappedMO->getName() << sendl;

    nbOutDOFs = mappedMO->getSize();
    nbInDOFs =  this->getMState()->getSize();

    mappedFFMatrix = new sofa::component::linearsolver::CompressedRowSparseMatrix<OutBlockType>;
    std::cout << "Size of output matrix: " << nbOutDOFs << std::endl;
    mappedFFMatrix->resizeBloc(nbOutDOFs,nbOutDOFs);

    mappedFFMatrixAccessor= new sofa::component::linearsolver::DefaultMultiMatrixAccessor;
    mappedFFMatrixAccessor->addMechanicalState(mappedMO);
    mappedFFMatrixAccessor->setGlobalMatrix(mappedFFMatrix);
    mappedFFMatrixAccessor->setupMatrices();

    timer = new CTime();

}

template<class DataTypes>
void MappedBeamToTetraForceField<DataTypes>::addForce(const core::MechanicalParams*  mparams, DataVecDeriv& , const DataVecCoord& , const DataVecDeriv& )
{    
    //std::cout << "AddForce in MappedFF" << std::endl;
    mparams->kFactor();
}

template<class DataTypes>
void MappedBeamToTetraForceField<DataTypes>::addDForce(const core::MechanicalParams*  mparams , DataVecDeriv& , const DataVecDeriv& )
{    
    mparams->kFactor();
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

    mappedFFMatrix->clear();
    mappedForceField->addKToMatrix(mparams,mappedFFMatrixAccessor);

    const sofa::component::linearsolver::CompressedRowSparseMatrix<MappingBlockType>* mappingMatrix;
    mappingMatrix = dynamic_cast<const sofa::component::linearsolver::CompressedRowSparseMatrix<MappingBlockType> *>(mapping->getJ());
    if (mappingMatrix == NULL) {
        serr << "Mapping " << mapping->getName() << " doesn't build matrix!" << sendl;
        return;
    }

    MappingBlockType mappingBuffer;
    OutBlockType outBlockBuffer;

    sofa::component::linearsolver::CompressedRowSparseMatrix<TempBlockType> tempMatrix;
    tempMatrix.resizeBloc(nbInDOFs, nbOutDOFs);

    //stopTime = (double)timer->getTime();
    //std::cout << "Length1 = " << stopTime - startTime << std::endl;
    //startTime = (double)timer->getTime();

    /*
     * K11 = mat
     * K22 = mappedFFMatrix
     *   J = mappingMatrix
     * */
    //operation K11 += Jt * K22 * J evaluated in two steps
    // NOTE: K11 is reffered to as IN and K22 as OUT matrix according to the
    //       direction of _mapping.

    //op1 : tempMatrix = Jt * K22
    for (int mappingRowIndex = 0; mappingRowIndex < mappingMatrix->nBlocRow; mappingRowIndex++) {   //through X, must take each row  (but test can be added)
        for (MappingColBlockConstIterator mappingColIter = mappingMatrix->bRowBegin(mappingRowIndex); mappingColIter < mappingMatrix->bRowEnd(mappingRowIndex); mappingColIter++) {  //take non zero blocks in row, determines the row in K)
            MappingBlockConstAccessor mappingBlock = mappingColIter.bloc();
            const MappingBlockType& mappingBlockData = *(const MappingBlockType*)mappingBlock.elements(mappingBuffer.ptr());
            int mappingColIndex = mappingBlock.getCol();
            for (OutColBlockConstIterator outColIter = mappedFFMatrix->bRowBegin(mappingRowIndex); outColIter < mappedFFMatrix->bRowEnd(mappingRowIndex); outColIter++) {
                OutBlockConstAccessor outBlock = outColIter.bloc();
                const OutBlockType& outBlockData = *(const OutBlockType*)outBlock.elements(outBlockBuffer.ptr());
                int outColIndex = outBlock.getCol();
                TempBlockType tempBlockData(0.0);
                //multiply the block, could be done more efficiently
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 6; j++)
                        for (int k = 0; k < 6; k++)
                            tempBlockData(i,j) += mappingBlockData(k,i)*outBlockData(k,j);
                tempMatrix.blocAdd(mappingColIndex, outColIndex, tempBlockData.ptr());
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
                InBlockType inBlockData(0.0);
                //multiply the block, could be done more efficiently
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++) {
                        for (int k = 0; k < 6; k++)
                            inBlockData(i,j) += tempBlockData(i,k) * mappingBlockData(k,j);
                        //because input matrix is  BaseMatrix, we do it like this for now
                        mat->add(3*tempRowIndex+i, 3*mappingColIndex+j, inBlockData(i,j));
                    }
                //mat->blocAdd(tempRowIndex,mappingColIndex,inBlockData.ptr());   //if mat is block matrix
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

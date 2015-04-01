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
 * RotationFinder.inl
 *
 *  Created on: 14 avr. 2009
 *      Author: froy
 */

#ifndef SOFA_COMPONENT_CONTAINER_SHAPEMATCHINGROTATIONFINDER_INL
#define SOFA_COMPONENT_CONTAINER_SHAPEMATCHINGROTATIONFINDER_INL

#include <sofa/component/container/ShapeMatchingRotationFinder.h>

#include <sofa/core/behavior/RotationMatrix.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/helper/decompose.h>
#include <sofa/helper/gl/template.h>


namespace sofa
{

namespace component
{

namespace container
{

using namespace sofa::core::topology;

template <class DataTypes>
ShapeMatchingRotationFinder<DataTypes>::ShapeMatchingRotationFinder()
: topo(NULL)
, axisToFlip(initData(&axisToFlip, int(-1), "axisToFlip", "Flip Axis"))
, showRotations(initData(&showRotations, bool(false), "showRotations", "Show Rotations"))
, neighborhoodLevel(initData(&neighborhoodLevel, int(1), "neighborhoodLevel", "Neighborhood level"))
, numOfClusters(initData(&numOfClusters, int(1), "numOfClusters", "Number of clusters"))
, maxIter(initData(&maxIter, unsigned(500), "maxIter", "Number of iterations to build the neighborhood"))
, epsilon(initData(&epsilon, Real(0.0000000001), "epsilon", "epsilon"))
, radius(initData(&radius, Real(0.001), "radius", "radius between Cm and point position"))
{

}

template <class DataTypes>
ShapeMatchingRotationFinder<DataTypes>::~ShapeMatchingRotationFinder()
{

}

template <class DataTypes>
void ShapeMatchingRotationFinder<DataTypes>::init()
{
        //Retrieve informations
	//- Mechanical State
	this->getContext()->get(mechanicalState);

	if (!mechanicalState)
	{
	    serr << "ERROR: RotationFinder has not found any MechanicalState" << sendl;
	    return;
	}

	//- Topology Container
 	this->getContext()->get(topo);
// 
// 	if (!topo)
// 	{
// 	    serr << "ERROR: RotationFinder has not found any BaseMeshTopology" << sendl;
// 	    return;
// 	}

	oldRestPositionSize = 0;
}

template <class DataTypes>
void ShapeMatchingRotationFinder<DataTypes>::ComputeNeighborhoodFromNeighborhood()
{
	Neighborhood neighborhood;
	for (unsigned int i=0; i<lastPointNeighborhood.size(); ++i)
	{
		neighborhood.clear();
		Neighborhood::const_iterator it, itEnd;
		for (it = lastPointNeighborhood[i].begin(), itEnd = lastPointNeighborhood[i].end(); it != itEnd ; ++it)
		{
			const helper::vector<Point>& vertexVertexShell = topo->getVerticesAroundVertex(*it);
			for (unsigned int j=0 ; j<vertexVertexShell.size(); ++j)
			{
				neighborhood.insert(vertexVertexShell[j]);
			}
		}

		std::insert_iterator<Neighborhood> result(pointNeighborhood[i], pointNeighborhood[i].begin());
		set_union(pointNeighborhood[i].begin(), pointNeighborhood[i].end(), neighborhood.begin(), neighborhood.end(), result);
		lastPointNeighborhood[i] = neighborhood;
	}
}

template <class DataTypes>
void ShapeMatchingRotationFinder<DataTypes>::computeNeighborhood()
{
    const unsigned int nbPoints =  mechanicalState->getSize();
    if(topo && neighborhoodLevel.getValue())
    {
        pointNeighborhood.resize(nbPoints);
        lastPointNeighborhood.resize(nbPoints);

        for (unsigned int i=0; i<nbPoints; ++i)
        {
            pointNeighborhood[i].insert(i);
            lastPointNeighborhood[i].insert(i);
        }

        for (int i=0; i<neighborhoodLevel.getValue() ; ++i)
        {
            ComputeNeighborhoodFromNeighborhood();
        }
    }
    else
    {
        const VecCoord& X0 = mechanicalState->read(sofa::core::ConstVecCoordId::restPosition())->getValue();
        helper::vector< unsigned int > clusterInPoint;
    
        clusterInPoint.resize(nbPoints);
        pointNeighborhood.resize(min<int>(numOfClusters.getValue(), nbPoints));
    
        if (pointNeighborhood.size()>1)
        {
            Coord cmObject;
    
            for(unsigned int i=0; i<nbPoints; ++i)
            {
                pointNeighborhood[i%pointNeighborhood.size()].insert(i);
                cmObject += X0[i];
            }
            cmObject /= nbPoints;
    
            Xcm0.clear();
            for(unsigned int i=0; i<pointNeighborhood.size(); ++i)
                Xcm0.push_back(cmObject);
            Xcm.resize(pointNeighborhood.size());
            
            unsigned int iter = 0;
            bool changed;
            do
            {
                changed = false;
                iter++;
                for(unsigned int i=0; i<pointNeighborhood.size(); ++i)
                {
                    Coord center;
                    Neighborhood::const_iterator it, itEnd;
                    for (it = pointNeighborhood[i].begin(), itEnd = pointNeighborhood[i].end(); it != itEnd ; ++it)
                    {
                        center += X0[*it];
                    }
                    if (pointNeighborhood[i].size())
                    {
                        center /= pointNeighborhood[i].size();
                        changed |= ((center - Xcm0[i]).norm2() > epsilon.getValue());
                        Xcm0[i] = center;
                        pointNeighborhood[i].clear();
                    }
                }
                double minDist, dist;
                int nearestCm0;
                for(unsigned int i=0; i<nbPoints; ++i)
                {
                    minDist = (Xcm0[0] - X0[i]).norm2();
                    nearestCm0 = 0;
    
                    for (unsigned int j=1; j<Xcm0.size(); ++j)
                    {
                        dist = (Xcm0[j] - X0[i]).norm2();
    
                        if (min<double>(dist, minDist) == dist)
                        {
                            minDist = dist;
                            nearestCm0 = j;
                        }
                    }
                    pointNeighborhood[nearestCm0].insert(i);
                    clusterInPoint[i] = nearestCm0;
                }
            }while(changed&&(iter<maxIter.getValue()));
            
            double minDist, dist;
            int nearestCm0;
            for(unsigned int i=0; i<nbPoints; ++i)
            {
                if (clusterInPoint[i] != 0)
                { 
                    minDist = (Xcm0[0] - X0[i]).norm2();
                    nearestCm0 = 0;
                }
                else
                {
                    minDist = (Xcm0[1] - X0[i]).norm2();
                    nearestCm0 = 1;
                }
    
                for (unsigned int j=nearestCm0+1; j<Xcm0.size(); ++j)
                {
                    if (clusterInPoint[i] != j)
                    { 
                        dist = (Xcm0[j] - X0[i]).norm2();
    
                        if (min<double>(dist, minDist) == dist)
                        {
                            minDist = dist;
                            nearestCm0 = j;
                        }
                    }
                }
    
                pointNeighborhood[nearestCm0].insert(i);
            }
            
            
            for(unsigned int i=0; i<nbPoints; ++i)
            {
                for(unsigned int j=0; j<pointNeighborhood.size(); ++j)
                {
                    if ((Xcm0[j] - X0[i]).norm2() < radius.getValue())
                    {
                        pointNeighborhood[j].insert(i);
                    }
                }
            }
        }
    }
}

template <class DataTypes>
void ShapeMatchingRotationFinder<DataTypes>::computeQT()
{

	Xcm.resize(pointNeighborhood.size());
	Xcm0.resize(pointNeighborhood.size());

        const VecCoord& restPositions = mechanicalState->read(sofa::core::ConstVecCoordId::restPosition())->getValue();

	for (unsigned int i=0;i<pointNeighborhood.size();++i)
	{
		Coord cm;
		Neighborhood::const_iterator it, itEnd;
		for (it = pointNeighborhood[i].begin(), itEnd = pointNeighborhood[i].end(); it != itEnd ; ++it)
			cm += restPositions[*it];
		cm /= pointNeighborhood[i].size();
		Xcm0[i] = cm;
	}
}

template <class DataTypes>
const typename ShapeMatchingRotationFinder<DataTypes>::VecNeighborhood& ShapeMatchingRotationFinder<DataTypes>::getNeighborhood()
{
    return pointNeighborhood;
}

template <class DataTypes>
void ShapeMatchingRotationFinder<DataTypes>::flipAxis(typename ShapeMatchingRotationFinder<DataTypes>::Mat3x3 & rotation)
{
	int axis = axisToFlip.getValue();
	if (axis >= 0 && axis <= 2)
	{
		for(unsigned int i=0 ; i<3;i++)
			rotation[i][axis] *= -1;
	}
}

template <class DataTypes>
const helper::vector<typename ShapeMatchingRotationFinder<DataTypes>::Mat3x3>& ShapeMatchingRotationFinder<DataTypes>::getRotations()
{
        const VecCoord& currentPositions = mechanicalState->read(sofa::core::ConstVecCoordId::position())->getValue();
        const VecCoord& restPositions = mechanicalState->read(sofa::core::ConstVecCoordId::restPosition())->getValue();

//	unsigned int nbPoints =  mechanicalState->getSize();

	if (currentPositions.size() < 3)
	{
                serr << "RotationFinder : problem with mechanical state; return ID matrix..." << sendl;
		rotations.clear();
		return rotations;
	}
	//if mechanical state has changed, we must compute again x0_cm and qT
	if(oldRestPositionSize != restPositions.size())
	{
		computeNeighborhood();
//                createClusters();
		computeQT();
		oldRestPositionSize = restPositions.size();
	}
	
        unsigned int nbShapes = pointNeighborhood.size();

	rotations.resize(nbShapes);
	dRotations.clear();
	dRotations_dx.clear();
	vecA.resize(nbShapes);

	for (unsigned int i=0 ; i<nbShapes ; i++)
	{
	    //we compute A_pq matrix
	    Mat3x3 A_pq(0);

	    Coord temp;
	    defaulttype::Mat<3,1,Real> p;
	    defaulttype::Mat<1,3,Real> qT;

            Coord center;
	    Neighborhood::const_iterator it, itEnd;
	    for (it = pointNeighborhood[i].begin(), itEnd = pointNeighborhood[i].end(); it != itEnd ; ++it)
		center += currentPositions[*it];

	    center /= pointNeighborhood[i].size();
	    Xcm[i] = center;

	    for (it = pointNeighborhood[i].begin(), itEnd = pointNeighborhood[i].end(); it != itEnd ; ++it)
	    {
		Coord neighbor = currentPositions[*it];
		
		temp = (neighbor - center);

		for (unsigned int k=0 ;k<3 ;k++)
			p[k][0] = temp[k];

		qT[0] = restPositions[*it] - Xcm0[i];

		A_pq += 1*(p*qT);
	    }
        Mat3x3 R;
        helper::Decompose<Real>::polarDecomposition(A_pq, R);

	    //test R is rotation or symmetry
	    //-> negative if symmetry
	    if (determinant(R) < 0)
		flipAxis(R);

	    vecA[i] = A_pq;
	    rotations[i] = R;
	}

	return rotations;
}

template <class DataTypes>
void ShapeMatchingRotationFinder<DataTypes>::getRotations(defaulttype::BaseMatrix * m,int offset) {
    if (component::linearsolver::RotationMatrix<Real> * diag = dynamic_cast<component::linearsolver::RotationMatrix<Real> *>(m)) {
        const VecCoord& currentPositions = mechanicalState->read(sofa::core::ConstVecCoordId::position())->getValue();
        const VecCoord& restPositions = mechanicalState->read(sofa::core::ConstVecCoordId::restPosition())->getValue();

//	unsigned int nbPoints =  mechanicalState->getSize();

	if (currentPositions.size() < 3)
	{
                serr << "RotationFinder : problem with mechanical state; return ID matrix..." << sendl;
		rotations.clear();
		return ;
	}
	//if mechanical state has changed, we must compute again x0_cm and qT
	if(oldRestPositionSize != restPositions.size())
	{
		computeNeighborhood();
//                createClusters();
		computeQT();
		oldRestPositionSize = restPositions.size();
	}
	
        unsigned int nbShapes = pointNeighborhood.size();

	diag->getVector().resize(nbShapes*9);

	for (unsigned int i=0 ; i<nbShapes ; i++)
	{
	    //we compute A_pq matrix
	    Mat3x3 A_pq(0);

	    Coord temp;
	    defaulttype::Mat<3,1,Real> p;
	    defaulttype::Mat<1,3,Real> qT;

            Coord center;
	    Neighborhood::const_iterator it, itEnd;
	    for (it = pointNeighborhood[i].begin(), itEnd = pointNeighborhood[i].end(); it != itEnd ; ++it)
		center += currentPositions[*it];

	    center /= pointNeighborhood[i].size();
	    Xcm[i] = center;

	    for (it = pointNeighborhood[i].begin(), itEnd = pointNeighborhood[i].end(); it != itEnd ; ++it)
	    {
		Coord neighbor = currentPositions[*it];
		
		temp = (neighbor - center);

		for (unsigned int k=0 ;k<3 ;k++)
			p[k][0] = temp[k];

		qT[0] = restPositions[*it] - Xcm0[i];

		A_pq += 1*(p*qT);
	    }
        Mat3x3 R;
        helper::Decompose<Real>::polarDecomposition(A_pq, R);

	    //test R is rotation or symmetry
	    //-> negative if symmetry
	    if (determinant(R) < 0)
		flipAxis(R);

	    *((Mat3x3 *) &diag->getVector()[i*9]) = R;
	}
    } else {
	getRotations();
	m->resize(rotations.size()*3,rotations.size()*3);
	
	for (unsigned i=0;i<rotations.size();i++) {
	    int e = i*3+offset;
	    m->set(e+0,e+0,rotations[i][0][0]);
	    m->set(e+0,e+1,rotations[i][0][1]);
	    m->set(e+0,e+2,rotations[i][0][2]);
	    
	    m->set(e+1,e+0,rotations[i][1][0]);
	    m->set(e+1,e+1,rotations[i][1][1]);
	    m->set(e+1,e+2,rotations[i][1][2]);
	    
	    m->set(e+2,e+0,rotations[i][2][0]);
	    m->set(e+2,e+1,rotations[i][2][1]);
	    m->set(e+2,e+2,rotations[i][2][2]);
	}
    }
}

template <class DataTypes>
const helper::vector<typename ShapeMatchingRotationFinder<DataTypes>::DMat3x3>& ShapeMatchingRotationFinder<DataTypes>::getDRotations()
{
	//const VecCoord& restPositions = *mechanicalState->getX0();
	//const VecDeriv& dx = *mechanicalState->getV();

	unsigned int nbShapes = pointNeighborhood.size();

	dRotations.resize(nbShapes);
	const Real epsilon = (Real) 0.000001;
	for (unsigned int i=0 ; i<nbShapes ; i++)
	{
	     DMat3x3 dR;
	    for (int l=0;l<3;++l)
	        for (int c=0;c<3;++c)
	    {
		Mat3x3 A = vecA[i];
		A[l][c]+=epsilon;
        Mat3x3 R;
        helper::Decompose<Real>::polarDecomposition(A, R);
	    //test R is rotation or symmetry
	    //-> negative if symmetry
	    if (determinant(R) < 0)
		flipAxis(R);
	     dR[l*3+c] = R;
	     dR[l*3+c] -= rotations[i];
	     dR[l*3+c] *= (((Real)1.0)/epsilon);
	    }
	    dRotations[i] = dR;
	}

	return dRotations;
}

template <class DataTypes>
const helper::vector<typename ShapeMatchingRotationFinder<DataTypes>::Mat3x3>& ShapeMatchingRotationFinder<DataTypes>::getDRotations(const VecDeriv& dx)
{
	if (dRotations.empty()) getDRotations();
        const VecCoord& restPositions = mechanicalState->read(sofa::core::ConstVecCoordId::restPosition())->getValue();

	unsigned int nbShapes = pointNeighborhood.size();

	dRotations_dx.resize(nbShapes);

	for (unsigned int i=0 ; i<nbShapes ; i++)
	{

	    //we compute dA_pq matrix
	    Mat3x3 dA_pq(0);

	    Coord temp;
	    defaulttype::Mat<3,1,Real> dp;
	    defaulttype::Mat<1,3,Real> qT;

            Coord dcenter;
	    Neighborhood::const_iterator it, itEnd;
	    for (it = pointNeighborhood[i].begin(), itEnd = pointNeighborhood[i].end(); it != itEnd ; ++it)
		dcenter += dx[*it];

	    dcenter /= pointNeighborhood[i].size();

	    for (it = pointNeighborhood[i].begin(), itEnd = pointNeighborhood[i].end(); it != itEnd ; ++it)
	    {
		Coord dneighbor = dx[*it];
		
		temp = (dneighbor - dcenter);

		for (unsigned int k=0 ;k<3 ;k++)
			dp[k][0] = temp[k];

		qT[0] = restPositions[*it] - Xcm0[i];

		dA_pq += 1*(dp*qT);
	    }
	    Mat3x3 dR;
	    for (int l=0;l<3;++l)
	        for (int c=0;c<3;++c)
		dR += dRotations[i][l*3+c] * dA_pq[l][c];

	    dRotations_dx[i] = dR;
	}
	return dRotations_dx;
}

template <class DataTypes>
void ShapeMatchingRotationFinder<DataTypes>::handleEvent(sofa::core::objectmodel::Event* event)
{
	 if (sofa::core::objectmodel::KeypressedEvent* ev = dynamic_cast<sofa::core::objectmodel::KeypressedEvent*>(event))
	    {
	        switch(ev->getKey())
	        {
				case 'i':
				case 'I':
					std::cout << "Rotation Matrix: " << std::endl;
					std::cout << getRotations() << std::endl;
					break;
			}
		}
}

template <class DataTypes>
void ShapeMatchingRotationFinder<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (showRotations.getValue())
    {
        const VecCoord& currentPositions = mechanicalState->read(sofa::core::ConstVecCoordId::position())->getValue();
	
        getRotations();

	glDisable(GL_LIGHTING);

	glBegin(GL_LINES);
	for (unsigned int i=0 ; i<rotations.size() ; i++)
	{
		glColor3f(1.0f,0.0f,0.0f);
		helper::gl::glVertexT(currentPositions[i]);
		helper::gl::glVertexT(currentPositions[i] + rotations[i].col(0));

		glColor3f(0.0f,1.0f,0.0f);
		helper::gl::glVertexT(currentPositions[i]);
		helper::gl::glVertexT(currentPositions[i] + rotations[i].col(1));

		glColor3f(0.0f,0.0f,1.0f);
		helper::gl::glVertexT(currentPositions[i]);
		helper::gl::glVertexT(currentPositions[i] + rotations[i].col(2));
	}
	glEnd();
        
        glEnable(GL_LIGHTING);
    }
        
    if (vparams->displayFlags().getShowForceFields())
    {
        const VecCoord& currentPositions = mechanicalState->read(sofa::core::ConstVecCoordId::position())->getValue();
        
        if (!showRotations.getValue())
            getRotations();
        
        glDisable(GL_LIGHTING);
        
        glBegin(GL_LINES);
        
        float r, g, b;
        
        for (unsigned int i=0 ; i<Xcm0.size() ; ++i)
        {
            r = (float)((i*7543)%11)/11;
            g = (float)((i*1357)%13)/13;
            b = (float)((i*4829)%17)/17;

            glColor3f(r,g,b);

            Neighborhood::const_iterator it, itEnd;
            for (it = pointNeighborhood[i].begin(), itEnd = pointNeighborhood[i].end(); it != itEnd ; ++it)
            {
                helper::gl::glVertexT(Xcm[i]);
                helper::gl::glVertexT(currentPositions[*it]);
            }
        }
        glEnd();

        glEnable(GL_LIGHTING);
    }
}

} // namespace container

} // namespace component

} // namespace sofa


#endif /* SOFA_COMPONENT_CONTAINER_SHAPEMATCHINGROTATIONFINDER_INL */

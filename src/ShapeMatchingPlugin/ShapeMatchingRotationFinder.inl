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

#include <ShapeMatchingPlugin/ShapeMatchingRotationFinder.h>

#include <sofa/linearalgebra/RotationMatrix.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/helper/decompose.h>
#include <ShapeMatchingPlugin/polar_decomposition_3x3.h>


namespace sofa::component::container
{

using namespace sofa::core::topology;

template <class DataTypes>
ShapeMatchingRotationFinder<DataTypes>::ShapeMatchingRotationFinder()
: l_mechanicalState(initLink("mechanicalState", "link to the mechanical state"))
, l_topology(initLink("topology", "link to the topology container"))
, d_axisToFlip(initData(&d_axisToFlip, int(-1), "axisToFlip", "Flip Axis"))
, d_showRotations(initData(&d_showRotations, bool(false), "showRotations", "Show Rotations"))
, d_neighborhoodLevel(initData(&d_neighborhoodLevel, int(1), "neighborhoodLevel", "Neighborhood level"))
, d_numOfClusters(initData(&d_numOfClusters, unsigned(1) , "numOfClusters", "Number of clusters"))
, d_maxIter(initData(&d_maxIter, unsigned(500), "maxIter", "Number of iterations to build the neighborhood"))
, d_epsilon(initData(&d_epsilon, Real(0.0000000001), "epsilon", "epsilon"))
, d_radius(initData(&d_radius, Real(0.001), "radius", "radius between Cm and point position"))
, m_oldRestPositionSize(0)
{

}

template <class DataTypes>
ShapeMatchingRotationFinder<DataTypes>::~ShapeMatchingRotationFinder()
{

}

template <class DataTypes>
void ShapeMatchingRotationFinder<DataTypes>::init()
{
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);

    //Retrieve informations
	//- Mechanical State
    if (!l_mechanicalState.get())
    {
        sofa::core::sptr<core::behavior::MechanicalState<DataTypes>> mechanicalState;
        this->getContext()->get(mechanicalState);
        if(!mechanicalState)
        {
            msg_error() << "did not found any Mechanical State.";
            return;
        }
        else
        {
            l_mechanicalState.set(mechanicalState);
        }
    }

	//- Topology Container
    if (!l_topology.get())
    {
        sofa::core::sptr<core::topology::BaseMeshTopology> topology;
        this->getContext()->get(topology);
        if (!topology)
        {
            msg_info() << "did not find any topology, will build standard point clusters.";
        }
        else
        {
            l_topology.set(topology);
        }
    }

    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}

template <class DataTypes>
void ShapeMatchingRotationFinder<DataTypes>::ComputeNeighborhoodFromNeighborhood()
{
	Neighborhood neighborhood;
	for (unsigned int i=0; i<m_lastPointNeighborhood.size(); ++i)
	{
		neighborhood.clear();
		Neighborhood::const_iterator it, itEnd;
		for (it = m_lastPointNeighborhood[i].begin(), itEnd = m_lastPointNeighborhood[i].end(); it != itEnd ; ++it)
		{
			const type::vector<Point>& vertexVertexShell = l_topology->getVerticesAroundVertex(*it);
			for (unsigned int j=0 ; j<vertexVertexShell.size(); ++j)
			{
				neighborhood.insert(vertexVertexShell[j]);
			}
		}

		std::insert_iterator<Neighborhood> result(m_pointNeighborhood[i], m_pointNeighborhood[i].begin());
		set_union(m_pointNeighborhood[i].begin(), m_pointNeighborhood[i].end(), neighborhood.begin(), neighborhood.end(), result);
		m_lastPointNeighborhood[i] = neighborhood;
	}
}

template <class DataTypes>
void ShapeMatchingRotationFinder<DataTypes>::computeNeighborhood()
{
    const unsigned int nbPoints =  l_mechanicalState->getSize();
    if(l_topology && d_neighborhoodLevel.getValue())
    {
        m_pointNeighborhood.resize(nbPoints);
        m_lastPointNeighborhood.resize(nbPoints);

        for (unsigned int i=0; i<nbPoints; ++i)
        {
            m_pointNeighborhood[i].insert(i);
            m_lastPointNeighborhood[i].insert(i);
        }

        for (int i=0; i<d_neighborhoodLevel.getValue() ; ++i)
        {
            ComputeNeighborhoodFromNeighborhood();
        }
    }
    else
    {
        const VecCoord& X0 = l_mechanicalState->read(sofa::core::vec_id::read_access::restPosition)->getValue();
        type::vector< unsigned int > clusterInPoint;
    
        clusterInPoint.resize(nbPoints);
        m_pointNeighborhood.resize(std::min(d_numOfClusters.getValue(), nbPoints));
    
        if (m_pointNeighborhood.size()>1)
        {
            Coord cmObject;
    
            for(unsigned int i=0; i<nbPoints; ++i)
            {
                m_pointNeighborhood[i%m_pointNeighborhood.size()].insert(i);
                cmObject += X0[i];
            }
            cmObject /= nbPoints;
    
            m_Xcm0.clear();
            for(unsigned int i=0; i<m_pointNeighborhood.size(); ++i)
                m_Xcm0.push_back(cmObject);
            m_Xcm.resize(m_pointNeighborhood.size());
            
            unsigned int iter = 0;
            bool changed;
            do
            {
                changed = false;
                iter++;
                for(unsigned int i=0; i<m_pointNeighborhood.size(); ++i)
                {
                    Coord center;
                    Neighborhood::const_iterator it, itEnd;
                    for (it = m_pointNeighborhood[i].begin(), itEnd = m_pointNeighborhood[i].end(); it != itEnd ; ++it)
                    {
                        center += X0[*it];
                    }
                    if (m_pointNeighborhood[i].size())
                    {
                        center /= m_pointNeighborhood[i].size();
                        changed |= ((center - m_Xcm0[i]).norm2() > d_epsilon.getValue());
                        m_Xcm0[i] = center;
                        m_pointNeighborhood[i].clear();
                    }
                }
                double minDist, dist;
                int nearestCm0;
                for(unsigned int i=0; i<nbPoints; ++i)
                {
                    minDist = (m_Xcm0[0] - X0[i]).norm2();
                    nearestCm0 = 0;
    
                    for (unsigned int j=1; j<m_Xcm0.size(); ++j)
                    {
                        dist = (m_Xcm0[j] - X0[i]).norm2();
    
                        if (std::min(dist, minDist) == dist)
                        {
                            minDist = dist;
                            nearestCm0 = j;
                        }
                    }
                    m_pointNeighborhood[nearestCm0].insert(i);
                    clusterInPoint[i] = nearestCm0;
                }
            }while(changed&&(iter<d_maxIter.getValue()));
            
            double minDist, dist;
            int nearestCm0;
            for(unsigned int i=0; i<nbPoints; ++i)
            {
                if (clusterInPoint[i] != 0)
                { 
                    minDist = (m_Xcm0[0] - X0[i]).norm2();
                    nearestCm0 = 0;
                }
                else
                {
                    minDist = (m_Xcm0[1] - X0[i]).norm2();
                    nearestCm0 = 1;
                }
    
                for (unsigned int j=nearestCm0+1; j<m_Xcm0.size(); ++j)
                {
                    if (clusterInPoint[i] != j)
                    { 
                        dist = (m_Xcm0[j] - X0[i]).norm2();
    
                        if (std::min(dist, minDist) == dist)
                        {
                            minDist = dist;
                            nearestCm0 = j;
                        }
                    }
                }
    
                m_pointNeighborhood[nearestCm0].insert(i);
            }
            
            
            for(unsigned int i=0; i<nbPoints; ++i)
            {
                for(unsigned int j=0; j<m_pointNeighborhood.size(); ++j)
                {
                    if ((m_Xcm0[j] - X0[i]).norm2() < d_radius.getValue())
                    {
                        m_pointNeighborhood[j].insert(i);
                    }
                }
            }
        }
    }
}

template <class DataTypes>
void ShapeMatchingRotationFinder<DataTypes>::computeQT()
{
	m_Xcm.resize(m_pointNeighborhood.size());
	m_Xcm0.resize(m_pointNeighborhood.size());

    const VecCoord& restPositions = l_mechanicalState->read(sofa::core::vec_id::read_access::restPosition)->getValue();

	for (unsigned int i=0;i<m_pointNeighborhood.size();++i)
	{
		Coord cm;
		Neighborhood::const_iterator it, itEnd;
		for (it = m_pointNeighborhood[i].begin(), itEnd = m_pointNeighborhood[i].end(); it != itEnd ; ++it)
			cm += restPositions[*it];
		cm /= m_pointNeighborhood[i].size();
		m_Xcm0[i] = cm;
	}
}

template <class DataTypes>
const typename ShapeMatchingRotationFinder<DataTypes>::VecNeighborhood& ShapeMatchingRotationFinder<DataTypes>::getNeighborhood()
{
    return m_pointNeighborhood;
}

template <class DataTypes>
void ShapeMatchingRotationFinder<DataTypes>::flipAxis(typename ShapeMatchingRotationFinder<DataTypes>::Mat3x3 & rotation)
{
	int axis = d_axisToFlip.getValue();
	if (axis >= 0 && axis <= 2)
	{
		for(unsigned int i=0 ; i<3;i++)
			rotation[i][axis] *= -1;
	}
}

template <class DataTypes>
const type::vector<typename ShapeMatchingRotationFinder<DataTypes>::Mat3x3>& ShapeMatchingRotationFinder<DataTypes>::getRotations()
{
    const VecCoord& currentPositions = l_mechanicalState->read(sofa::core::vec_id::read_access::position)->getValue();
    const VecCoord& restPositions = l_mechanicalState->read(sofa::core::vec_id::read_access::restPosition)->getValue();

	if (currentPositions.size() < 3)
	{
		msg_error() << "RotationFinder : problem with mechanical state; return ID matrix...";
		m_rotations.clear();
		return m_rotations;
	}
	//if mechanical state has changed, we must compute again x0_cm and qT
	if(m_oldRestPositionSize != restPositions.size())
	{
		computeNeighborhood();
		computeQT();
		m_oldRestPositionSize = restPositions.size();
	}
	
    auto nbShapes = m_pointNeighborhood.size();

	m_rotations.resize(nbShapes);
	m_dRotations.clear();
	m_dRotations_dx.clear();
	m_vecA.resize(nbShapes);

	for (unsigned int i=0 ; i<nbShapes ; i++)
	{
	    //we compute A_pq matrix
	    Mat3x3 A_pq(0);

	    Coord temp;
	    type::Mat<3,1,Real> p;
	    type::Mat<1,3,Real> qT;

            Coord center;
	    Neighborhood::const_iterator it, itEnd;
	    for (it = m_pointNeighborhood[i].begin(), itEnd = m_pointNeighborhood[i].end(); it != itEnd ; ++it)
		center += currentPositions[*it];

	    center /= m_pointNeighborhood[i].size();
	    m_Xcm[i] = center;

	    for (it = m_pointNeighborhood[i].begin(), itEnd = m_pointNeighborhood[i].end(); it != itEnd ; ++it)
	    {
		Coord neighbor = currentPositions[*it];
		
		temp = (neighbor - center);

		for (unsigned int k=0 ;k<3 ;k++)
			p[k][0] = temp[k];

		qT[0] = restPositions[*it] - m_Xcm0[i];

		A_pq += 1*(p*qT);
	    }
        Mat3x3 R, H;
		polar::polar_decomposition(R.ptr(), H.ptr(), A_pq.ptr());

	    //test R is rotation or symmetry
	    //-> negative if symmetry
	    if (determinant(R) < 0)
		flipAxis(R);

	    m_vecA[i] = A_pq;
	    m_rotations[i] = R;
	}

	return m_rotations;
}

template <class DataTypes>
void ShapeMatchingRotationFinder<DataTypes>::getRotations(linearalgebra::BaseMatrix * m,int offset) 
{
    if (auto* diag = dynamic_cast<linearalgebra::RotationMatrix<Real> *>(m))
    {
        const VecCoord& currentPositions = l_mechanicalState->read(sofa::core::vec_id::read_access::position)->getValue();
        const VecCoord& restPositions = l_mechanicalState->read(sofa::core::vec_id::read_access::restPosition)->getValue();

        if (currentPositions.size() < 3)
        {
            msg_error() << "RotationFinder : problem with mechanical state; return ID matrix...";
            m_rotations.clear();
            return;
        }

        //if mechanical state has changed, we must compute again x0_cm and qT
        if (m_oldRestPositionSize != restPositions.size())
        {
            computeNeighborhood();
            computeQT();
            m_oldRestPositionSize = restPositions.size();
        }

        const auto nbShapes = m_pointNeighborhood.size();

        diag->getVector().resize(nbShapes * 9);

        for (unsigned int i = 0; i < nbShapes; i++)
        {
            //we compute A_pq matrix
            Mat3x3 A_pq(0);

            Coord temp;
            type::Mat<3, 1, Real> p;
            type::Mat<1, 3, Real> qT;

            Coord center;
            Neighborhood::const_iterator it, itEnd;
            for (it = m_pointNeighborhood[i].begin(), itEnd = m_pointNeighborhood[i].end(); it != itEnd; ++it)
                center += currentPositions[*it];

            center /= m_pointNeighborhood[i].size();
            m_Xcm[i] = center;

            for (it = m_pointNeighborhood[i].begin(), itEnd = m_pointNeighborhood[i].end(); it != itEnd; ++it)
            {
                Coord neighbor = currentPositions[*it];

                temp = (neighbor - center);

                for (unsigned int k = 0; k < 3; k++)
                    p[k][0] = temp[k];

                qT[0] = restPositions[*it] - m_Xcm0[i];

                A_pq += 1 * (p * qT);
            }
            Mat3x3 R;
            helper::Decompose<Real>::polarDecomposition(A_pq, R);

            //test R is rotation or symmetry
            //-> negative if symmetry
            if (determinant(R) < 0)
                flipAxis(R);

            *((Mat3x3*)&diag->getVector()[i * 9]) = R;
        }
    }
    else 
    {
        getRotations();
        m->resize(sofa::SignedIndex(m_rotations.size() * 3), sofa::SignedIndex(m_rotations.size() * 3));

        for (unsigned i = 0; i < m_rotations.size(); i++) {
            int e = i * 3 + offset;
            m->set(e + 0, e + 0, m_rotations[i][0][0]);
            m->set(e + 0, e + 1, m_rotations[i][0][1]);
            m->set(e + 0, e + 2, m_rotations[i][0][2]);

            m->set(e + 1, e + 0, m_rotations[i][1][0]);
            m->set(e + 1, e + 1, m_rotations[i][1][1]);
            m->set(e + 1, e + 2, m_rotations[i][1][2]);

            m->set(e + 2, e + 0, m_rotations[i][2][0]);
            m->set(e + 2, e + 1, m_rotations[i][2][1]);
            m->set(e + 2, e + 2, m_rotations[i][2][2]);
        }
    }
}

template <class DataTypes>
const type::vector<typename ShapeMatchingRotationFinder<DataTypes>::DMat3x3>& ShapeMatchingRotationFinder<DataTypes>::getDRotations()
{
	auto nbShapes = m_pointNeighborhood.size();

	m_dRotations.resize(nbShapes);
	constexpr Real epsilon = (Real) 0.000001;
	constexpr auto invEpsilon = (((Real)1.0) / epsilon);
	for (unsigned int i=0 ; i<nbShapes ; i++)
	{
	    for (int l=0;l<3;++l)
	        for (int c=0;c<3;++c)
	    {
			Mat3x3 A = m_vecA[i];
			A[l][c]+=epsilon;
			Mat3x3 R, H;
			polar::polar_decomposition(R.ptr(), H.ptr(), A.ptr());
			//test R is rotation or symmetry
			//-> negative if symmetry
			if (determinant(R) < 0)
				flipAxis(R);
			m_dRotations[i][l * 3 + c] = R;
			m_dRotations[i][l * 3 + c] -= m_rotations[i];
			m_dRotations[i][l * 3 + c] *= invEpsilon;
	    }
	}

	return m_dRotations;
}

template <class DataTypes>
const type::vector<typename ShapeMatchingRotationFinder<DataTypes>::Mat3x3>& ShapeMatchingRotationFinder<DataTypes>::getDRotations(const VecDeriv& dx)
{
	if (m_dRotations.empty())
	{
		getDRotations();
	}

    const VecCoord& restPositions = l_mechanicalState->read(sofa::core::vec_id::read_access::restPosition)->getValue();

	const auto nbShapes = m_pointNeighborhood.size();

	m_dRotations_dx.resize(nbShapes);

	Mat3x3 dA_pq(0);

	Coord temp;
	type::Mat<3, 1, Real> dp;
	type::Mat<1, 3, Real> qT;
	Coord dcenter;

	for (unsigned int i=0 ; i<nbShapes ; i++)
	{
	    //we compute dA_pq matrix
		dA_pq.clear();
		dcenter.clear();

		for (const auto& pn : m_pointNeighborhood[i])
		{
			dcenter += dx[pn];
		}

	    dcenter /= m_pointNeighborhood[i].size();

		for (const auto& pn : m_pointNeighborhood[i])
	    {
			const Coord& dneighbor = dx[pn];
			temp = (dneighbor - dcenter);

			for (unsigned int k = 0; k < 3; k++)
			{
				dp[k][0] = temp[k];
			}

			qT[0] = restPositions[pn] - m_Xcm0[i];
			dA_pq += (dp*qT);
	    }
	    
		m_dRotations_dx[i].clear();
		for (int l = 0; l < 3; ++l)
		{
			for (int c = 0; c < 3; ++c)
			{
				m_dRotations_dx[i] += m_dRotations[i][l * 3 + c] * dA_pq[l][c];
			}
		}
	}
	return m_dRotations_dx;
}

template <class DataTypes>
void ShapeMatchingRotationFinder<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
	vparams->drawTool()->saveLastState();

    if (d_showRotations.getValue())
    {

        const VecCoord& currentPositions = l_mechanicalState->read(sofa::core::vec_id::read_access::position)->getValue();
	
        getRotations();

		vparams->drawTool()->disableLighting();

		std::vector<type::Vec3> vertices;
		std::vector<type::RGBAColor> colors;

		for (unsigned int i=0 ; i<m_rotations.size() ; i++)
		{
			vertices.push_back(currentPositions[i]);
			vertices.push_back(currentPositions[i] + m_rotations[i].col(0));
			colors.push_back(type::RGBAColor::red());

			vertices.push_back(currentPositions[i]);
			vertices.push_back(currentPositions[i] + m_rotations[i].col(1));
			colors.push_back(type::RGBAColor::green());

			vertices.push_back(currentPositions[i]);
			vertices.push_back(currentPositions[i] + m_rotations[i].col(2));
			colors.push_back(type::RGBAColor::blue());
		}

		vparams->drawTool()->drawLines(vertices, 1.0f, colors);

		vparams->drawTool()->enableLighting();
    }
        
    if (vparams->displayFlags().getShowForceFields())
    {
        const VecCoord& currentPositions = l_mechanicalState->read(sofa::core::vec_id::read_access::position)->getValue();
        
        if (!d_showRotations.getValue())
            getRotations();
        
		vparams->drawTool()->disableLighting();

		std::vector<type::Vec3> vertices;
		std::vector<type::RGBAColor> colors;
        
        float r, g, b;
        
        for (unsigned int i=0 ; i<m_Xcm0.size() ; ++i)
        {
            r = (float)((i*7543)%11)/11;
            g = (float)((i*1357)%13)/13;
            b = (float)((i*4829)%17)/17;

            Neighborhood::const_iterator it, itEnd;
            for (it = m_pointNeighborhood[i].cbegin(), itEnd = m_pointNeighborhood[i].cend(); it != itEnd ; ++it)
            {
				vertices.push_back(m_Xcm[i]);
				vertices.push_back(currentPositions[*it]);

				colors.push_back({ r, g, b, 1.0f });
            }
        }

		vparams->drawTool()->drawLines(vertices, 1.0f, colors);

		vparams->drawTool()->enableLighting();
    }

	vparams->drawTool()->restoreLastState();
}

} // namespace sofa::component::container

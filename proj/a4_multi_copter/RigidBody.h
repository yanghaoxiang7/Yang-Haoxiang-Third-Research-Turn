#ifndef __RigidBody_h__
#define __RigidBody_h__
#include "Common.h"

//////////////////////////////////////////////////////////////////////////
////Rigid body simulator
template<int d> class RigidBody
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using MatrixD=Matrix<real,d>;
public:
	VectorD position=VectorD::Zero();
	VectorD velocity=VectorD::Zero();
	MatrixD R=MatrixD::Identity();
	MatrixD Rt=MatrixD::Identity();
	VectorD omega=VectorD::Zero();

	const VectorD WorldVectorToLocalVector(const VectorD& v){return Rt * v;}
	const VectorD LocalVectorToWorldVector(const VectorD& v){return R * v;}
	const VectorD WorldPointToLocalPoint(const VectorD& p){return Rt * (p - position);}
	const VectorD LocalPointToWorldPoint(const VectorD& p){return position + R * p;}
};

#endif
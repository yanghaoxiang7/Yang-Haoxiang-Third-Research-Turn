//#####################################################################
// Mass-spring deformable model
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################

#ifndef __SoftBodyMassSpring_h__
#define __SoftBodyMassSpring_h__
#include "Common.h"
#include "Particles.h"

template<int d> class SoftBodyMassSpring
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using MatrixD=Matrix<real,d>;
public:
	////Spring parameters
	Particles<d> particles;
	Array<Vector2i> springs;
	Array<real> rest_length;
	Array<real> ks;
	Array<real> kd;

	////Boundary nodes
	Hashtable<int,VectorD> boundary_nodes;

	////Body force
	VectorD g=VectorD::Unit(1)*(real)-1.;
	
	enum class TimeIntegration{ExplicitEuler,ImplicitEuler} time_integration=TimeIntegration::ExplicitEuler;

	////Implicit time integration
	SparseMatrixT K;
	VectorX u,b;

	virtual void Initialize()
	{
		////Initialize default spring parameters for standard tests
		real ks_0=(real)1,kd_0=(real)1;
		switch(time_integration){
		case TimeIntegration::ExplicitEuler:{
			ks_0=(real)1e2;
			kd_0=(real)1e1;
		}break;
		case TimeIntegration::ImplicitEuler:{
			ks_0=(real)1e5;
			kd_0=(real)1e1;			
		}break;}

		////Allocate arrays for springs and parameters
		rest_length.resize(springs.size());
		for(int i=0;i<(int)springs.size();i++){const Vector2i& s=springs[i];
			rest_length[i]=(particles.X(s[0])-particles.X(s[1])).norm();}
		ks.resize(springs.size(),ks_0);
		kd.resize(springs.size(),kd_0);

		////Allocate sparse matrix if using implicit time integration 
		////This function needs to be called for only once since the mesh doesn't change during the simulation)
		if(time_integration==TimeIntegration::ImplicitEuler)
			Initialize_Implicit_K_And_b();
	}

	virtual void Advance(const real dt)
	{
		switch(time_integration){
		case TimeIntegration::ExplicitEuler:
			Advance_Explicit_Euler(dt);break;
		case TimeIntegration::ImplicitEuler:
			Advance_Implicit_Euler(dt);break;}
	}
	
	////Set boundary nodes
	void Set_Boundary_Node(const int p,const VectorD v=VectorD::Zero()){boundary_nodes[p]=v;}
	
	bool Is_Boundary_Node(const int p){return boundary_nodes.find(p)!=boundary_nodes.end();}
	
	void Enforce_Boundary_Conditions()
	{
		for(auto p:boundary_nodes){
			int idx=p.first;					////get boundary particle index
			const VectorD& v=p.second;			////get boundary particle velocity
			particles.V(idx)=v;					////set boundary particle velocity
			particles.F(idx)=VectorD::Zero();}	////clear boundary particle force
	}

	//////////////////////////////////////////////////////////////////////////
	////P1 TASK: explicit Euler integration and spring force calculation

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P1 TASK): explicit Euler time integration 
	void Advance_Explicit_Euler(const real dt)
	{
		Particle_Force_Accumulation();
		////Update particle velocity and position
		for (int i = 0; i < (int)particles.Size(); i++)
		{
			VectorD x = particles.X(i), v = particles.V(i);
			particles.V(i) += dt / particles.M(i) * particles.F(i);
			particles.X(i) += dt * v;
		}
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P1 TASK): compute spring force f_ij=f_s+f_d 
	VectorD Spring_Force_Calculation(const int idx)
	{
		int i = springs[idx][0], j = springs[idx][1];
		real k_s = ks[idx], k_d = kd[idx], l_0 = rest_length[idx];
		VectorD x_i = particles.X(i), x_j = particles.X(j);
		VectorD v_i = particles.V(i), v_j = particles.V(j);
		
		real nm = (x_j - x_i).norm();
		VectorD f_s = k_s * (nm - l_0) * ((x_j - x_i) / nm);
		VectorD f_d = k_d * ((v_j - v_i).dot((x_j - x_i) / nm)) * ((x_j - x_i) / nm);
//		printf("%lf %lf %lf\n", (double) nm, (double) f_s.norm(), (double) f_d.norm());
		
		VectorD f_ij = f_s + f_d;
		return f_ij;
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P1 TASK): accumulate spring forces to particles
	void Particle_Force_Accumulation()
	{
		////Clear forces on particles
		for(int i=0;i<particles.Size();i++){particles.F(i)=VectorD::Zero();}

		////Accumulate body forces
		for(int i=0;i<particles.Size();i++){
			particles.F(i)+=particles.M(i)*g;}

		////Accumulate spring forces
		for (int idx = 0; idx < (int)springs.size(); idx++)
		{
			VectorD f = Spring_Force_Calculation(idx);
			int i = springs[idx][0], j = springs[idx][1];
			particles.F(i) += f;
			particles.F(j) -= f;
		}

		////Enforce boundary conditions
		Enforce_Boundary_Conditions();
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P2 TASK): 
	////Construct K, step 1: initialize the matrix structure 
	void Initialize_Implicit_K_And_b()
	{
		int n=d*particles.Size();
		K.resize(n,n);u.resize(n);u.fill((real)0);b.resize(n);b.fill((real)0);
		Array<TripletT> elements;
		for(int s=0;s<(int)springs.size();s++){int i=springs[s][0];int j=springs[s][1];
			Add_Block_Triplet_Helper(i,i,elements);
			Add_Block_Triplet_Helper(i,j,elements);
			Add_Block_Triplet_Helper(j,i,elements);
			Add_Block_Triplet_Helper(j,j,elements);}
		K.setFromTriplets(elements.begin(),elements.end());
		K.makeCompressed();	
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P2 TASK): 
	////Construct K, step 2: fill nonzero elements in K
	void Update_Implicit_K_And_b(const real dt)
	{
		////Clear K and b
		K.setZero();
		b.fill((real)0);

		int n = particles.Size(), m = (int) springs.size();

		// M
		for (int i = 0; i < n; i++)
			SparseFunc::Add_Block<d,MatrixD>(K,i,i,particles.M(i) * MatrixD::Identity());
		
		// -dt * D (in K) and -dt * D * v_n (in b)
		MatrixD Kd;
		for (int s = 0; s < m; s++)
		{
			Compute_Kd_Block(s, Kd);
			int i = springs[s][0], j = springs[s][1];
			Add_Block_Helper(K, i, j, Kd * (-dt));
			
			VectorD tmp;
			tmp = Kd * particles.V(i) + (-Kd) * particles.V(j);
			Add_Block(b, i, tmp * (-dt));
			tmp = (-Kd) * particles.V(i) + Kd * particles.V(j);
			Add_Block(b, j, tmp * (-dt));
		}
		
		// -dt^2 * K
		MatrixD Ks;
		for (int s = 0; s < m; s++)
		{
			Compute_Ks_Block(s, Ks);
			int i = springs[s][0], j = springs[s][1];
			Add_Block_Helper(K, i, j, Ks * (-dt * dt));
		}
		
		// Mv
		for (int i = 0; i < n; i++)
			Add_Block(b, i, particles.M(i) * particles.V(i));
		
		// dt * f
		for (int i = 0; i < n; i++)
			Add_Block(b, i, dt * particles.F(i));
	}

	//////////////////////////////////////////////////////////////////////////
	////P2 TASK: Implicit Euler time integration

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P2 TASK): 
	////Construct K, step 2.1: compute spring force derivative
	void Compute_Ks_Block(const int s, MatrixD& Ks)
	{
		int i = springs[s][0], j = springs[s][1];
		VectorD d = particles.X(j) - particles.X(i);
		real l_0 = rest_length[s];
		Ks = (l_0 / d.norm() - (real) 1) * MatrixD::Identity();
		Ks -= l_0 / (d.norm() * d.norm() * d.norm())
			* d * d.transpose();
		Ks = Ks * ks[s];
	}

	//////////////////////////////////////////////////////////////////////////
	////YOUR IMPLEMENTATION (P2 TASK): 
	////Construct K, step 2.2: compute damping force derivative
	void Compute_Kd_Block(const int s, MatrixD& Kd)
	{
		int i = springs[s][0], j = springs[s][1];
		VectorD d = particles.X(j) - particles.X(i);
		Kd = d * d.transpose() / (d.norm() * d.norm());
		Kd = -Kd * kd[s];
	}

	////Implicit Euler time integration
	void Advance_Implicit_Euler(const real dt)
	{
		Particle_Force_Accumulation();
		Update_Implicit_K_And_b(dt);

		for(int i=0;i<particles.Size();i++){
			for(int j=0;j<d;j++)u[i*d+j]=particles.V(i)[j];}	////set initial guess to be the velocity from the last time step

		SparseSolver::CG(K,u,b);	////solve Ku=b using Conjugate Gradient

		for(int i=0;i<particles.Size();i++){
			VectorD v;for(int j=0;j<d;j++)v[j]=u[i*d+j];
			particles.V(i)=v;
			particles.X(i)+=particles.V(i)*dt;
		}
	}

protected:
	////Add block nonzeros to sparse matrix elements (for initialization)
	void Add_Block_Triplet_Helper(const int i,const int j,Array<TripletT>& elements)
	{for(int ii=0;ii<d;ii++)for(int jj=0;jj<d;jj++)elements.push_back(TripletT(i*d+ii,j*d+jj,(real)0));}

	////Add block Ks to K_ij
	void Add_Block_Helper(SparseMatrixT& K,const int i,const int j,const MatrixD& Ks)
	{
		SparseFunc::Add_Block<d,MatrixD>(K,i,i,Ks);
		SparseFunc::Add_Block<d,MatrixD>(K,j,j,Ks);
		if(!Is_Boundary_Node(i)&&!Is_Boundary_Node(j)){
			SparseFunc::Add_Block<d,MatrixD>(K,i,j,-Ks);
			SparseFunc::Add_Block<d,MatrixD>(K,j,i,-Ks);}
	}

	////Set block values on a vector
	void Set_Block(VectorX& b,const int i,const VectorD& bi)
	{for(int ii=0;ii<d;ii++)b[i*d+ii]=bi[ii];}

	////Add block values to a vector
	void Add_Block(VectorX& b,const int i,const VectorD& bi)
	{for(int ii=0;ii<d;ii++)b[i*d+ii]+=bi[ii];}
};

#endif

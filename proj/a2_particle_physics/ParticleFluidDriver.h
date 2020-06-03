//#####################################################################
// Particle Fluid Driver
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################
#ifndef __ParticleFluidDriver_h__
#define __ParticleFluidDriver_h__
#include <random>
#include "Common.h"
#include "Driver.h"
#include "OpenGLMarkerObjects.h"
#include "OpenGLCommon.h"
#include "OpenGLWindow.h"
#include "OpenGLViewer.h"
#include "ParticleFluid.h"

#define Times 0.5

template<int d> class ParticleFluidDriver : public Driver, public OpenGLViewer
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using Base=Driver;
	real dt=.02;
	ParticleFluid<d> fluid;
	Array<OpenGLSolidCircle*> opengl_circles;

	Bowl<d>* bowl=nullptr;
public:
	virtual void Initialize()
	{
		////driver initialization, initialize simulation data
		real dx=.35 * Times;int nx=int(20 / Times);int ny=int(20 / Times);
		for(int i=0;i<nx;i++){for(int j=0;j<ny;j++){
			VectorD pos;pos[0]=(real)i*dx-1.;pos[1]=(real)j*dx+3.;
			Add_Particle(pos);}}

		bowl=new Bowl<d>(VectorD::Unit(1)*8,8);
		fluid.env_objects.push_back(bowl);

		fluid.Initialize();

		////viewer initialization, initialize visualization data
		OpenGLViewer::Initialize();
	}

	////synchronize simulation data to visualization data, called in OpenGLViewer::Initialize()
	virtual void Initialize_Data()
	{
		if(bowl){
			auto opengl_circle=Add_Interactive_Object<OpenGLCircle>();
			opengl_circle->n=64;
			opengl_circle->pos=V3(bowl->center);
			opengl_circle->radius=bowl->radius;
			opengl_circle->color=OpenGLColor(1.f,.6f,.2f);
			opengl_circle->line_width=4.f;
			opengl_circle->Set_Data_Refreshed();
			opengl_circle->Initialize();}

		for(int i=0;i<fluid.particles.Size();i++){
			Add_Solid_Circle(i);}
	}

	void Sync_Simulation_And_Visualization_Data()
	{
		for(int i=0;i<fluid.particles.Size();i++){
			auto opengl_circle=opengl_circles[i];
			opengl_circle->pos=V3(fluid.particles.X(i));
			opengl_circle->Set_Data_Refreshed();}
	}

	////update simulation and visualization for each time step
	virtual void Toggle_Next_Frame()
	{
		fluid.Advance(dt);
		Sync_Simulation_And_Visualization_Data();
		OpenGLViewer::Toggle_Next_Frame();
	}

	virtual void Run()
	{
		OpenGLViewer::Run();
	}

	////User interaction
	virtual bool Mouse_Click(int left,int right,int mid,int x,int y,int w,int h)
	{
		if(left!=1){return false;}
		Vector3f win_pos=opengl_window->Project(Vector3f::Zero());
		Vector3f pos=opengl_window->Unproject(Vector3f((float)x,(float)y,win_pos[2]));
		VectorD p_pos;for(int i=0;i<d;i++)p_pos[i]=(real)pos[i];
		real r=.1*static_cast<float>(rand()%1000)/1000.+.15;
		Add_Particle(p_pos);
		Add_Solid_Circle(fluid.particles.Size()-1);
		return true;
	}

protected:
	void Add_Particle(VectorD pos,real m=1.)
	{
		int i=fluid.particles.Add_Element();	////return the last element's index
		fluid.particles.X(i)=pos;
		fluid.particles.V(i)=VectorD::Zero();
		fluid.particles.R(i)=0.3 * Times; // 0.1
		fluid.particles.M(i)=m;
		fluid.particles.D(i)=1. * Times;
	}

	void Add_Solid_Circle(const int i)
	{
		OpenGLColor c(0.5f,0.5f,0.5f,1.f);
		auto opengl_circle=Add_Interactive_Object<OpenGLSolidCircle>();
		opengl_circles.push_back(opengl_circle);
		opengl_circle->pos=V3(fluid.particles.X(i));
		opengl_circle->radius=fluid.particles.R(i);
		opengl_circle->color=c;
		opengl_circle->Set_Data_Refreshed();
		opengl_circle->Initialize();	
	}


	////Helper function to convert a vector to 3d, for c++ template
	Vector3 V3(const Vector2& v2){return Vector3(v2[0],v2[1],.0);}
	Vector3 V3(const Vector3& v3){return v3;}
};
#endif

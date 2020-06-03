//#####################################################################
// Particle Sand Driver
// Dartmouth COSC 89.18/189.02: Computational Methods for Physical Systems, Assignment starter code
// Contact: Bo Zhu (bo.zhu@dartmouth.edu)
//#####################################################################
#ifndef __ParticleSandDriver_h__
#define __ParticleSandDriver_h__
#include <random>
#include "Common.h"
#include "Driver.h"
#include "OpenGLMarkerObjects.h"
#include "OpenGLCommon.h"
#include "OpenGLWindow.h"
#include "OpenGLViewer.h"
#include "ParticleSand.h"

template<int d> class ParticleSandDriver : public Driver, public OpenGLViewer
{using VectorD=Vector<real,d>;using VectorDi=Vector<int,d>;using Base=Driver;
	real dt=.02;
	ParticleSand<d> sand;
	Array<OpenGLSolidCircle*> opengl_circles;

	Bowl<d>* bowl=nullptr;
public:
	virtual void Initialize()
	{
		////driver initialization, initialize simulation data
		Add_Particle(VectorD::Unit(1)*1.,.4,1.);
		Add_Particle(VectorD::Unit(1)*2.,.2,1.);
		Add_Particle(VectorD::Unit(0)*.2+VectorD::Unit(1)*2.5,.1,1.);
		Add_Particle(VectorD::Unit(0)*-.2+VectorD::Unit(1)*2.5,.3,1.);

		bowl=new Bowl<d>(VectorD::Unit(1)*8,8);
		sand.env_objects.push_back(bowl);

		////viewer initialization, initialize visualization data
		OpenGLViewer::Initialize();
	}

	////synchronize simulation data to visualization data
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

		for(int i=0;i<sand.particles.Size();i++){
			Add_Solid_Circle(i);}
	}

	void Sync_Simulation_And_Visualization_Data()
	{
		for(int i=0;i<sand.particles.Size();i++){
			auto opengl_circle=opengl_circles[i];
			opengl_circle->pos=V3(sand.particles.X(i));
			opengl_circle->Set_Data_Refreshed();}
	}

	////update simulation and visualization for each time step
	virtual void Toggle_Next_Frame()
	{
		sand.Advance(dt);
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
		Add_Particle(p_pos,r,1.);
		Add_Solid_Circle(sand.particles.Size()-1);
		return true;
	}

protected:
	void Add_Particle(VectorD pos,real r=.1,real m=1.)
	{
		int i=sand.particles.Add_Element();	////return the last element's index
		sand.particles.X(i)=pos;
		sand.particles.V(i)=VectorD::Zero();
		sand.particles.R(i)=r;
		sand.particles.M(i)=m;
	}

	void Add_Solid_Circle(const int i)
	{
		OpenGLColor c;
		for(int i=0;i<3;i++){
			c.rgba[i]=static_cast<float>(rand()%1000)/1000.f;}
		auto opengl_circle=Add_Interactive_Object<OpenGLSolidCircle>();
		opengl_circles.push_back(opengl_circle);
		opengl_circle->pos=V3(sand.particles.X(i));
		opengl_circle->radius=sand.particles.R(i);
		opengl_circle->color=c;
		opengl_circle->Set_Data_Refreshed();
		opengl_circle->Initialize();	
	}

	////Helper function to convert a vector to 3d, for c++ template
	Vector3 V3(const Vector2& v2){return Vector3(v2[0],v2[1],.0);}
	Vector3 V3(const Vector3& v3){return v3;}
};
#endif

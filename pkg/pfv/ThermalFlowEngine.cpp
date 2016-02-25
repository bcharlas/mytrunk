 
/*************************************************************************
*  Copyright (C) 2014 by Bruno Chareyre <bruno.chareyre@hmg.inpg.fr>     *
*                                                                        *
*  This program is free software; it is licensed under the terms of the  *
*  GNU General Public License v2 or later. See file LICENSE for details. *
*************************************************************************/

// This is an example of how to derive a new FlowEngine with additional data and possibly completely new behaviour.
// Every functions of the base engine can be overloaded, and new functions can be added

//keep this #ifdef as long as you don't really want to realize a final version publicly, it will save compilation time for everyone else
//when you want it compiled, you can pass -DDUMMYFLOW to cmake, or just uncomment the following line
// #define DUMMYFLOW
#ifdef THERMALFLOW

#include "FlowEngine_ThermalFlowEngineT.hpp"
#include "FlowEngine.hpp"

/// We can add data to the Info types by inheritance
class ThermalCellInfo : public FlowCellInfo_ThermalFlowEngineT
{
	public:
	//
	double molarMass;
	molarMassGas=17.05;//in g/mol
	
	//private:
	double cellTemperature;
	cellTemperature=273;
	double cellCv;
	double cellCp;
	//double cellGamma;
	
	
	
	
};

class ThermalVertexInfo : public ThermalVertexInfo_ThermalFlowEngineT {
	public:
	
	//
	bool isThermal;
	isThermal=false;
	double bodyTemperature;
	bodyTemperature=273;
	double incTemperature;// increment of temperature
	incTemperature=0.;
	//double bodyCp;
	//double bodyHeatCond;// heat conductivity of the body
	//
	// sorption properties
	bool isNH3;
	isNH3 = false;
	double molarMass;
	molarMass=158.53;//in g/mol
	double enthalpySorption;
	enthalpySorption = 41400.0 // enthalpy in J/mol
	double entropySorption;
	entropySorption = 228.8 // entropy in J/mol/K
	double coefSorption;
	coefSorption =9.5*1e3 // mol/s
	double activationNRGSorption;
	activationNRGSorption=32.9* 1e3 // J/mol
	double specificR;
	specificR=8.31/molarMass //J/g/K
	double refPressure;
	refPressure=1.0// reference pressure in Pa
	double absorbedMass;
	absorbedMass=0.0//
	//double absorbedMassMax;//max mass of ammonia absorbable by the body
	double incMass;// increment of absorbed mass
	incMass=0.0;//
	double diffusionInSolid;//diffusion coefficient of NH3 between salt particles
	diffusionInSolid = 1e-3;
	//
	//void getBodyVolume();
	//void setBodyVolume();
	
	//private:
	//double bodyVolume;
	
};

typedef TemplateFlowEngine_ThermalFlowEngineT<ThermalCellInfo,ThermalVertexInfo> ThermalFlowEngineT;
REGISTER_SERIALIZABLE(ThermalFlowEngineT);
YADE_PLUGIN((ThermallowEngineT));

class ThermalFlowEngine : public ThermalFlowEngineT
{
	public :
	//We can overload every functions of the base engine to make it behave differently
	//if we overload action() like this, this engine is doing nothing in a standard timestep, it can still have useful functions
	virtual void action() {};
	
	//If a new function is specific to the derived engine, put it here, else go to the base TemplateFlowEngine
	//if it is useful for everyone
	void fancyFunction(Real what);

	YADE_CLASS_BASE_DOC_ATTRS_INIT_CTOR_PY(ThermalFlowEngine,ThermalFlowEngineT,"documentation here",
	,/*DummyFlowEngineT()*/,
	,
	// cell properties
	/*
	.def("setCellTemperature",&ThermalFlowEngine::setCellTemperature,(boost::python::arg("id"),boost::python::arg("cellTemperature"),"set temperature properties in cell 'id'.")
	.def("getCellTemperature",&ThermalFlowEngine::getCellTemperature,(boost::python::arg("id")),"get thermal properties in cell 'id'.")
	//
	.def("setCellCv",&ThermalFlowEngine::setCellCv,(boost::python::arg("id"),boost::python::arg("cellCv"),"set thermal properties in cell 'id'.")
	.def("getCellCv",&ThermalFlowEngine::getCellCv,(boost::python::arg("id")),"get thermal properties in cell 'id'.")
	.def("setCellCp",&ThermalFlowEngine::setCellCp,(boost::python::arg("id"),boost::python::arg("cellCp"),"set thermal properties in cell 'id'.")
	.def("getCellCp",&ThermalFlowEngine::getCellCp,(boost::python::arg("id")),"get thermal properties in cell 'id'.")
	//
	.def("initCellProperties",&ThermalFlowEngine::initCellProperties,(boost::python::arg("cellTemperature"), boost::python::arg("cellCv"), boost::python::arg("cellCp"),"set thermal properties in cells 'id'.")
	//
	//
	// Bodies properties
	//.def("setBodyTemperature",&ThermalFlowEngine::setBodyTemperature,(boost::python::arg("id"),boost::python::arg("bodyTemperature"),"set temperature properties of the body 'id'.")
	//.def("getBodyTemperature",&ThermalFlowEngine::getBodyTemperature,(boost::python::arg("id")),"get thermal properties in cell 'id'.")
	//
	
	.def("initBodiesProperties",&ThermalFlowEngine::initBodiesProperties,(boost::python::arg("listBodies"), 
									      boost::python::arg("bodyTemperature")
									      ,"set initial thermal properties in vertices.")
	
	*/
	//.def("setBodyVolume",&ThermalFlowEngine::setCellTemperature,(boost::python::arg("id"),boost::python::arg("bodyVolume")),"set volume of body 'id'.")
	//.def("getBodyTemperature",&ThermalFlowEngine::getCellTemperature,(boost::python::arg("id")),"get temperature in cell 'id'.")
	)
	DECLARE_LOGGER;
};
REGISTER_SERIALIZABLE(ThermalFlowEngine);

/*
void ThermalFlowEngine::setCellTemperature(int id, double temperature)
{
	solver->T[solver->currentTes].cellHandles[id]->info().cellTemperature = temperature;
}

double ThermalFlowEngine::getCellTemperature(int id)
{
	return solver->T[solver->currentTes].cellHandles[id]->info().cellTemperature;
}
//
void ThermalFlowEngine::setCellCv(int id, double Cv)
{
	solver->T[solver->currentTes].cellHandles[id]->info().cellCv = Cv;
}

double ThermalFlowEngine::setCellCv(int id)
{
	return solver->T[solver->currentTes].cellHandles[id]->info().cellCv;
}
//
void ThermalFlowEngine::setCellCp(int id, double Cp)
{
	solver->T[solver->currentTes].cellHandles[id]->info().cellCp = Cp;
}

double ThermalFlowEngine::setCellCp(int id)
{
	return solver->T[solver->currentTes].cellHandles[id]->info().cellCp;
}
//
void ThermalFlowEngine::initCellProperties(double temperature, double Cv, double Cp)
{
	FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
	{
	 cell->info().cellTemperature = temperature;
	 cell->info().cellCv = Cv;
	 cell->info().cellCp = Cp;
	}
}

//


void ThermalFlowEngine::setBodyTemperature(int id, double temperature)
{
	solver->T[solver->currentTes].cellHandles->vertex(i)->info().cellTemperature = temperature;
}

double ThermalFlowEngine::getBodyTemperature(int id)
{
	return solver->T[solver->currentTes].cellHandles[id]->info().cellTemperature;
}

void ThermalFlowEngine::initBodiesProperties(boost::python::list listBodies,  double temperature)
{
	FOREACH(CellHandle& cell, solver->T[solver->currentTes].cellHandles)
	{
	  for (unsigned int i=0;i<4;i++)
	  {
	    if cell->vertex(i)->info().id() boost::python::in listBodies
	    {
	      cell->vertex(i)->info().isThermal=true;
	      cell->vertex(i)->info().isNH3=true;
	      //
	      cell->vertex(i)->info().bodyTemperature=temperature;
	      
	      T[currentTes].vertexHandles[vhi.id()]->info().forces // IMPORTANT
	    }
	      
	    cell->vertex(i)->info().id()
	    
	    
	  }
	}
}
*/

solver->T[solver->currentTes].cellHandles[id]->vertex(i)->info().id()

YADE_PLUGIN((ThermalFlowEngine));

#endif //ThermalFLOW
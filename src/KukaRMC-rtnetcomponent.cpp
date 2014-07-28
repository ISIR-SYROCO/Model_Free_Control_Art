// Filename:  kukaMingXingController-rtnetcomponent.cpp
// Copyright: 2014 ISIR-CNRS
// Author:  Sovan Hak (hak@isir.upmc.fr) 
// Description:  

#include "KukaRMC-rtnetcomponent.hpp"
#include <rtt/Component.hpp>
#include <iostream>
#include <fstream>
#include <Eigen/LU>
#include <sstream>

#define PI	3.141592653589793238462643383280 /* pi */

KukaRMControllerRTNET::KukaRMControllerRTNET(std::string const& name) : m_dt(0.001), m_time(0.0), m_Fidx(0), m_FBoucle(false), m_delta(30), m_FSat(200.0), m_SL(30), m_L(10), m_FL(2), m_fastLen(1), m_slowLen(1), m_writeInFile(false), m_compteur(0),  FriRTNetExampleAbstract(name)
{
    this->addOperation("setFreq", &KukaRMControllerRTNET::setFreq, this, RTT::OwnThread);
    this->addOperation("setAmp", &KukaRMControllerRTNET::setAmp, this, RTT::OwnThread);
    this->addOperation("setBias", &KukaRMControllerRTNET::setBias, this, RTT::OwnThread);
    this->addOperation("setName", &KukaRMControllerRTNET::setName, this, RTT::OwnThread);
    this->addOperation("setAllowWritingData", &KukaRMControllerRTNET::setAllowWritingData, this, RTT::OwnThread);
    this->addOperation("setMAFDelta", &KukaRMControllerRTNET::setMAFDelta, this, RTT::OwnThread);
    this->addOperation("setAMAFLength", &KukaRMControllerRTNET::setAMAFLength, this, RTT::OwnThread);
    this->addOperation("setJointImpedance", &KukaRMControllerRTNET::setJointImpedance, this, RTT::OwnThread);
    
	/* Resizing all vectors */
	m_Kp.resize(LWRDOF);
    m_Kd.resize(LWRDOF);
    m_pid.resize(LWRDOF);
    m_q_des.resize(LWRDOF);
	m_qp_des.resize(LWRDOF);
	m_qpp_des.resize(LWRDOF);
	m_q_mes.resize(LWRDOF);
	m_qp_mes.resize(LWRDOF);
	m_qpp_mes.resize(LWRDOF);
	m_tau_mes.resize(LWRDOF);
	m_tau_FRI.resize(LWRDOF);
	m_amp.resize(LWRDOF);
	m_freq.resize(LWRDOF);
	m_bias.resize(LWRDOF);
	m_H.resize(LWRDOF,LWRDOF);
	m_c.resize(LWRDOF);
	m_g.resize(LWRDOF);
	m_FTab.resize(LWRDOF,m_delta);
	m_bufSum.resize(LWRDOF);
	m_prec.resize(LWRDOF);
	m_AMA.resize(LWRDOF);
	m_buf.resize(LWRDOF,m_L);

	/* Gains du PID (articulaire) */
	m_Kp << 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0;
	m_Kd << 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7;

	/* Initialisation a zero */
	m_pid.setZero();
	m_q_des.setZero();
	m_qp_des.setZero();
	m_qpp_des.setZero();
	m_qp_mes.setZero();
	m_qpp_mes.setZero();
	m_tau_mes.setZero();
	m_tau_FRI.setZero();
	m_bufSum.setZero();
	m_FTab.setZero();
	m_prec.setZero();
	m_AMA.setZero();
	m_buf.setZero();
    
    /* Modele calcule */
    m_model = new kukafixed("kuka");
    
    /* Dynamique du modele */
    m_H = m_model->getInertiaMatrix();
    m_c = m_model->getNonLinearTerms();
    m_g = m_model->getGravityTerms();
    
    /* Output Files */
    m_fluxErrName = "";
    m_fluxTrqName = "";
}

KukaRMControllerRTNET::~KukaRMControllerRTNET(){
	if(m_writeInFile)
	{
		m_fluxErr.seekp(-4,std::ios::end);
		m_fluxErr << "];%";
		m_fluxErr.close();
		m_fluxTrq.seekp(-4,std::ios::end);
		m_fluxTrq << "];%";
		m_fluxTrq.close();
		m_flux.close();
	}
	delete m_model;
	m_model = 0;
}

void KukaRMControllerRTNET::updateHook()
{
    std::string fri_mode("e_fri_unkown_mode");
    bool fri_cmd_mode = false;
    RTT::FlowStatus fs_event = iport_events.read(fri_mode);
    if (fri_mode == "e_fri_cmd_mode")
        fri_cmd_mode = true;
    else if (fri_mode == "e_fri_mon_mode")
        fri_cmd_mode = false;

    Eigen::VectorXd tau_dyn(LWRDOF);
    Eigen::VectorXd F(LWRDOF);
	
	std::vector<double> jntPos(LWRDOF);
	std::vector<double> jntVel(LWRDOF);
	std::vector<double> jntTrq(LWRDOF);
	std::vector<double> tau_FRI(LWRDOF);
	
	/* Get des mesures des capteurs */
	RTT::FlowStatus joint_pos_fs = iport_msr_joint_pos.read(jntPos);
	RTT::FlowStatus joint_vel_fs = iport_msr_joint_vel.read(jntVel);
	RTT::FlowStatus joint_trq_fs = iport_msr_joint_trq.read(jntTrq);
	
	if(fri_cmd_mode){
		if(m_compteur<5000)
		{
		/* Mise a jour des mesures */
		if(joint_trq_fs==RTT::NewData){
			for(unsigned int i=0;i<LWRDOF;i++){
				m_tau_mes[i] = jntTrq[i];
			}
		}
		if(joint_pos_fs==RTT::NewData){
			for(unsigned int i=0;i<LWRDOF;i++){
				m_q_mes[i] = jntPos[i];
			}
			m_model->setJointPositions(m_q_mes); //Mise a jour des positions dans le modele
		}
		if(joint_vel_fs==RTT::NewData){
			fdm(jntVel); //calcule l acceleration qpp
			for(unsigned int i=0;i<LWRDOF;i++){
				m_qp_mes[i] = jntVel[i];
			}
			m_model->setJointVelocities(m_qp_mes); //Mise a jour des vitesses dans le modele
		}
	
		/* Mise a jour des attributs */
		tau_dyn = m_tau_mes-m_tau_FRI;
		m_H = m_model->getInertiaMatrix();
		m_c = m_model->getNonLinearTerms();
		m_g = m_model->getGravityTerms();
		
		/* Estimation de F */
		F = m_qpp_mes-m_H.inverse()*(m_tau_FRI-m_c-m_g-m_pid);
		//F = m_qpp_mes-m_H.inverse()*(m_tau_FRI+tau_dyn-m_c-m_g-m_pid);
		estimateF(F);
		//estimateAMAF(F);
		
		/* Calcul des positions, vitesses et accelerations desirees */
		m_time = m_time+m_dt; //Mise a jour du temps
		nextDesVal();
		
		/* Calcul du PID */
		pidCalc();
		
		/* Calcul de la commande */
		m_tau_FRI.setZero(); //m_H*(m_qpp_des-F)+m_c+m_g+m_pid;
		//m_tau_FRI = m_H*(m_qpp_des-F)+m_c+m_g-tau_dyn+m_pid;
	
		if(m_writeInFile)
		{
			m_flux << std::endl;
			m_flux << "temps : " << m_time << std::endl;
			m_flux << "q_des : " << m_q_des.transpose() << std::endl;
			m_flux << "q_mes : " << m_q_mes.transpose() << std::endl;
			m_flux << "qp_des : " << m_qp_des.transpose() << std::endl;
			m_flux << "qp_mes : " << m_qp_mes.transpose() << std::endl;
			m_flux << "qpp_des : " << m_qpp_des.transpose() << std::endl;
			m_flux << "qpp_mes : " << m_qpp_mes.transpose() << std::endl;
			//m_flux << "F_est_sat : " << F.transpose() << std::endl;
			//m_flux << "PID : " << m_pid.transpose() << std::endl;
			//m_flux << "Hqpp : " << (m_H*m_qpp_des).transpose() << std::endl;
			//m_flux << "coriolis : " << m_c.transpose() << std::endl;
			//m_flux << "gravity : " << m_g.transpose() << std::endl;
			m_flux << "tau_FRI : " << m_tau_FRI.transpose() << std::endl;
			m_flux << "f_dynamics : " << tau_dyn.transpose() << std::endl;
			m_flux << "gravity : " << m_g.transpose() << std::endl;
			//m_flux << "tau_mes : " << m_tau_mes.transpose() << std::endl;
			m_fluxTrq << std::endl << m_tau_mes.transpose() << ";...";
			m_fluxErr << std::endl << (m_q_des-m_q_mes).transpose() << ";...";
		}
		
		for(unsigned int i=0;i<LWRDOF;i++){
			tau_FRI[i] = m_tau_FRI[i];
		}
		oport_add_joint_trq.write(tau_FRI);
		oport_joint_position.write(jntPos);
		}
		else{
			m_compteur++;
		}
	}
}


void KukaRMControllerRTNET::nextDesVal()
{
	Eigen::VectorXd period = m_freq*m_time;
	Eigen::VectorXd amp = Eigen::VectorXd::Zero(LWRDOF);
	
	if(m_time<1.0){
		amp = m_amp*m_time;
	}
	else{
		amp = m_amp;
	}
	for(unsigned int i=0;i<LWRDOF;i++)
	{
		m_q_des(i) = amp(i)*(period.array().sin())(i)+m_bias(i);
		m_qp_des(i) = amp(i)*m_freq(i)*(period.array().cos())(i);
		m_qpp_des(i) = -amp(i)*m_freq(i)*m_freq(i)*(period.array().sin())(i);
	}
}

void KukaRMControllerRTNET::estimateF(Eigen::VectorXd& F)
{
	for(unsigned int i=0;i<LWRDOF;i++)
	{
		m_bufSum(i) = m_bufSum(i)-m_FTab(i,m_Fidx);
		m_FTab(i,m_Fidx) = F(i);
		m_bufSum(i) = m_bufSum(i)+F(i);
	}
	
	m_Fidx = m_Fidx+1;
	if(m_Fidx>m_delta-1){
		m_Fidx = 0;
		m_FBoucle = true;
	}
	
	if(m_FBoucle){
		F = m_bufSum/m_delta;
	}
	else{
		F = m_bufSum/m_Fidx;
	}
	
	for(unsigned int i=0;i<LWRDOF;i++)
	{
		if((F.array().abs())(i)>m_FSat)
		{
			if(F(i)>0.0){
				F(i) = m_FSat;
			}
			else{
				F(i) = -m_FSat;
			}
		}
	}
}

void KukaRMControllerRTNET::estimateAMAF(Eigen::VectorXd& F)
{
	Eigen::VectorXd er(LWRDOF);
	Eigen::VectorXd sc(LWRDOF);
	Eigen::VectorXd precL(LWRDOF);
	double fast = 2/(m_fastLen+1);
	double slow = 2/(m_slowLen+1);
	
	for(unsigned int i=0;i<LWRDOF;i++)
	{
		m_bufSum(i) = m_bufSum(i)-m_FTab(i,m_Fidx);
		
		/* Memorisation de la nbAxe-m_lg valeur */
    	if(m_FBoucle){
        	precL(i) = m_buf(i,m_Fidx);
        }
    	else{
        	precL(i) = m_buf(i,0);
    	}
    
    	/* Mis a jour de la i-ieme valeur */
    	m_buf(i,m_Fidx) = ((F-m_prec).array().abs())(i);
    	
    	/* Mise a jour de la somme */
    	m_bufSum(i) = m_bufSum(i)+m_buf(i,m_Fidx); 
    	
    	/* Calcul de ER */
    	if(m_bufSum(i)==0){
        	er(i) = 0;
        }
    	else{
        	er(i) = ((F-precL).array().abs())(i)/m_bufSum(i);
    	}
    	
    	/* Calcul de sc */
    	sc(i) = (er(i)*(fast-slow)+slow)*(er(i)*(fast-slow)+slow);
    
    	/* Calcul de la moyenne mobile */
    	m_AMA(i) = m_AMA(i)*(1.0-sc(i))+sc(i)*F(i);
    	
    	/* Sortie de la fonction */
    	F(i) = m_AMA(i);
    	
    	m_prec(i) = m_buf(i,m_Fidx);
	}
	
	/* Mise a jour des parametres pour la prochaine utilisation */
    m_Fidx=m_Fidx+1;
    if(m_Fidx>m_L)
    {
        m_Fidx = 1;
        m_FBoucle = 1;
    }
    if(m_fastLen<m_FL){
        m_fastLen = m_fastLen+1;
    }
    if(m_slowLen<m_SL){
        m_slowLen = m_slowLen+1;
    }
}

void KukaRMControllerRTNET::fdm(std::vector<double> const& qp_new)
{
	for(unsigned int i=0;i<LWRDOF;i++){
		m_qpp_mes[i] = (qp_new[i]-m_qp_mes[i])/m_dt;
	}
}

void KukaRMControllerRTNET::pidCalc()
{
	for(unsigned int i=0;i<LWRDOF;i++){
		m_pid(i) = m_Kp(i)*(m_q_des(i)-m_q_mes(i))+m_Kd(i)*(m_qp_des(i)-m_qp_mes(i));
	}
}

void KukaRMControllerRTNET::setFreq(std::vector<double> const& freq)
{
	for(unsigned int i=0;i<LWRDOF;i++){
		m_freq[i] = freq[i];
	}
}

void KukaRMControllerRTNET::setAmp(std::vector<double> const& amp)
{
	for(unsigned int i=0;i<LWRDOF;i++){
		m_amp[i] = amp[i];
	}
}
void KukaRMControllerRTNET::setBias(std::vector<double> const& bias)
{
	for(unsigned int i=0;i<LWRDOF;i++){
		m_bias[i] = bias[i];
		m_q_mes[i] = m_bias[i];
	}
}

void KukaRMControllerRTNET::setName(std::string const& fluxErrName, std::string const& fluxTrqName)
{
	m_fluxErrName = fluxErrName;
	m_fluxTrqName = fluxTrqName;
}

void KukaRMControllerRTNET::setAllowWritingData(bool writeInFile)
{
	std::string ss;	
	
	m_writeInFile = writeInFile;
	if(m_writeInFile)
	{
		if(!(m_fluxErrName==""))
		{
			ss = "/home/kuka/src/groovy_workspace/model_free_control_art/"+m_fluxErrName;
			m_fluxErr.open(ss.c_str());
    		m_fluxErr << "err_q = ["; 
    		ss = "/home/kuka/src/groovy_workspace/model_free_control_art/"+m_fluxTrqName;
    		m_fluxTrq.open(ss.c_str());
    		m_fluxTrq << "tau_mes = [..."; 
    		m_flux.open("/home/kuka/src/groovy_workspace/model_free_control_art/data.dat");
    	}
    	else{
    		std::cout << "Error ! Name of files has not been provided" << std::endl;
    		return;
    	}
    }
}

void KukaRMControllerRTNET::setMAFDelta(unsigned const int delta)
{
	m_delta = delta;
	m_FTab.resize(LWRDOF,m_delta);
	m_FTab.setZero();
}

void KukaRMControllerRTNET::setAMAFLength(unsigned int const& FL, unsigned int const& L, unsigned int const& SL)
{
	m_SL = SL;
	m_FL = FL;
	m_L  =  L;
	m_buf.resize(LWRDOF,m_L);
	m_buf.setZero();
}

void KukaRMControllerRTNET::setJointImpedance(std::vector<double> &stiffness, std::vector<double> &damping)
{
	if(stiffness.size() != LWRDOF || damping.size() != LWRDOF)
	{
		std::cout << "Wrong vector size, should be " << LWRDOF << ", " << LWRDOF << std::endl;
		return;
	}
	else{
		lwr_fri::FriJointImpedance joint_impedance_command;
		for(unsigned int i = 0; i < LWRDOF; i++){
			joint_impedance_command.stiffness[i] = stiffness[i];
			joint_impedance_command.damping[i] = damping[i];
		}
		oport_joint_impedance.write(joint_impedance_command);
	}
}

/*
 * Using this macro, only one component may live
 * in one library *and* you may *not* link this library
 * with another component library. Use
 * ORO_CREATE_COMPONENT_TYPE()
 * ORO_LIST_COMPONENT_TYPE(kuka_ming_xing_controller)
 * In case you want to link with another library that
 * already contains components.
 *
 * If you have put your component class
 * in a namespace, don't forget to add it here too:
 */
ORO_CREATE_COMPONENT(KukaRMControllerRTNET)

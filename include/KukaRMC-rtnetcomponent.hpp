// Filename:  KukaRMC-rtnetcomponent.hpp
// Copyright: 2014 ISIR-CNRS
// Author:  Vincent SAY 
// Description: Orocos component using RTNET to compute the kuka model
//              from fri data and using a restricted model controller.

#ifndef KUKA_KUKA_RMC_RTNET_COMPONENT_HPP
#define KUKA_KUKA_RMC_RTNET_COMPONENT_HPP

#include <friRTNetExampleAbstract.hpp>
#include "kukafixed.h"
#include <Eigen/Dense>
#include <fstream>

class KukaRMControllerRTNET : public FriRTNetExampleAbstract{
	private:
		/* modele */
		kukafixed *m_model;
		
		/* Vecteurs des coordonnees */
		Eigen::VectorXd m_q_mes;
		Eigen::VectorXd m_qp_mes;
		Eigen::VectorXd m_qpp_mes;
		Eigen::VectorXd m_q_des;
		Eigen::VectorXd m_qp_des;
		Eigen::VectorXd m_qpp_des;
		
		/* Vecteurs des efforts */
		Eigen::VectorXd m_tau_mes;
		Eigen::VectorXd m_tau_FRI;
		Eigen::VectorXd m_tau_max;		
		
		/* Vecteurs des PID */
		Eigen::VectorXd m_pid;
		Eigen::VectorXd m_Kp;
		Eigen::VectorXd m_Kd;
		
		/* Vecteurs des signaux desires */
		Eigen::VectorXd m_amp;
		Eigen::VectorXd m_freq;
		Eigen::VectorXd m_bias;
		
		/* Matrice et vecteurs dynamiques */
		Eigen::MatrixXd m_H;
		Eigen::VectorXd m_c;
		Eigen::VectorXd m_g;
		
		/* Variables pour estimation de F */
		Eigen::MatrixXd m_FTab;
		Eigen::VectorXd m_bufSum;
		double m_FSat;
		unsigned int m_Fidx;
		bool m_FBoucle;
		unsigned int m_delta; //Fenetre de moyennage
		/* Variables pour AMAF */
		unsigned int m_SL;
		unsigned int m_L;
		unsigned int m_FL;
		Eigen::MatrixXd m_buf;
		Eigen::VectorXd m_prec;
		Eigen::VectorXd m_AMA;
		double m_fastLen;
		double m_slowLen;
		
		double m_dt; //Temps echantillonnage
		double m_time; //Temps actuel
		bool m_writeInFile;
		unsigned int m_compteur;
		std::ofstream m_fluxErr;
		std::string m_fluxErrName;
		std::ofstream m_fluxTrq;
		std::string m_fluxTrqName;
		std::ofstream m_flux;
		
	public:
		KukaRMControllerRTNET(std::string const& name);
		~KukaRMControllerRTNET();
		void updateHook();
		void nextDesVal();
		void estimateF(Eigen::VectorXd& F);
		void estimateAMAF(Eigen::VectorXd& F);
		void fdm(std::vector<double> const& qp_new);
		void pidCalc();
		
		void setFreq(std::vector<double> const& freq);
		void setAmp(std::vector<double> const& amp);
		void setBias(std::vector<double> const& bias);
		void setName(std::string const& fluxErrName, std::string const& fluxTrqName);
		void setAllowWritingData(bool writeInFile);
		void setMAFDelta(unsigned const int delta);
		void setAMAFLength(unsigned int const& FL, unsigned int const& L, unsigned int const& SL);
		void setJointImpedance(std::vector<double> &stiffness, std::vector<double> &damping);
};

#endif

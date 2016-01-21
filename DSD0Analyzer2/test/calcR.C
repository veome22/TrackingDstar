void calcR(float dataKpi, float dataKpiErr, float dataK3pi, float dataK3piErr, float MCKpi, float MCK3pi) {

int nMCK3piTot = 2425449;
int nMCKpiTot = 2225449; // FIXME

effKpi = MCKpi / nMCKpiTot;
effK3pi = MCK3pi / nMCK3piTot;

float R = (dataK3pi/dataKpi)*(effKpi/effK3pi);
float errR = sqrt( (dataKpiErr/dataKpi)**2 + (dataK3piErr/dataK3pi)**2 + (1./MCKpi) + (1./MCK3pi) + (1./nMCK3piTot) + (1./nMCKpiTot) );

std::cout<<"R = "<<R<<" +/- "<<errR<<std::endl;

float RPDG = 2.08;
float RPDGErr = 0.05;

effRatio = sqrt(R/RPDG);
effRatioErr = 0.5 * sqrt( (errR/R)**2 + (RPDGErr/RPDG)**2 );

std::cout<<"eff(Data)/eff(MC) = "<<effRatio<<" +/- "<<effRatioErr<<std::endl;
}

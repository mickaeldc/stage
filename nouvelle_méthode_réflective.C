#include <cmath>
#include <iostream>
#include <TMath.h>
#include <cmath>
#include "TF2.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"
#include "TGraph.h"
using namespace std;




void DataGenerator(double* RunPositionX,double* RunPositionY,int Nruns,double r)
{

	 TF2 *Gauss2D = new TF2("Gauss2D","gaus(x,[0],[1],[2])*gaus(y,[3],[4],[5])",-10,10,-10,10);
    	 TF2 *Source = new TF2("Source","gaus(x,[0],[6],[7])*gaus(y,[3],[8],[9])*gaus(x,[0],[1],[2])*gaus(y,[3],[4],[5])",-10,10,-10,10);
	 Gauss2D->SetParameter(0,1.);
   	 Gauss2D->SetParameter(2,1.353);
   	 Gauss2D->SetParameter(3,1.);
         Gauss2D->SetParameter(5,1.353);
         Source-> SetParameter(6,0);
         Source-> SetParameter(7,r);    	 
	 Source-> SetParameter(8,0);	
	 Source-> SetParameter(9,r);
	 Source->SetParameter(0,1.);	
   	 Source->SetParameter(3,1.);
 	 Source->SetParameter(2,1.353);
         Source->SetParameter(5,1.353);
	 
	 double PosXEvent = 0;
    	 double PosYEvent = 0;
    	 double EventNumber = 0;
    	 double RunNumber = 0;
    	 //double indice = 0;
    	TFile* fileTree = new TFile("test.root","RECREATE");
    	TTree* Events = new TTree("DSTevents","DSTevents");
    
    	Events->Branch("PosXEvent",&PosXEvent,"PosXEvent/D");
    	Events->Branch("PosYEvent",&PosYEvent,"PosYEvent/D");
    	Events->Branch("EventNumber",&EventNumber,"EventNumber/D");
    	Events->Branch("RunNumber",&RunNumber,"RunNumber/D");
	//Events->Branch("indice",&indice,"indice/D");
	
	int FrequenceBruit = 265;    	
	int TimePerRun =  28*60;	
	double NEventPerRunBruit = (FrequenceBruit * TimePerRun)/10;
    	double NEventSource = 0;
	int EventID = 0;
    
    	for (int j = 0; j < Nruns; j++)
	{
		Gauss2D->SetParameter(1,RunPositionX[j]);
	        Gauss2D->SetParameter(4,RunPositionY[j]);       	 	
		Source->SetParameter(1,RunPositionX[j]);
	        Source->SetParameter(4,RunPositionY[j]); 		
		TString foncBruit ="FoncBruit";	
		foncBruit += j+1;
		TString fonc = "FoncSource";
		fonc += j+1;	
		//TCanvas* FoncSource = new TCanvas(fonc,fonc);
		//Source -> Draw("surf1");		
		//TCanvas* Foncruit = new TCanvas(foncBruit,foncBruit);		
		//Gauss2D ->Draw("surf1");		
		RunNumber = j;
        	//cout << " Simulate Run number " << j+1 << " at position " << RunPositionX[j] << " " << RunPositionY[j] << endl;       
		for (int i = 0; i < (NEventPerRunBruit+NEventSource); i++)
		{
	    		if (i < NEventPerRunBruit)
			{
				//cout<<"Evenement de bruit numeros: "<< i<<endl;				
				EventNumber = i;
            			Gauss2D->GetRandom2(PosXEvent,PosYEvent);
				//indice = 0;
            			Events->Fill();
            			EventID++;
				                
			}
		
			if (i > NEventPerRunBruit)
			{
				//cout<<"Evenement de source numeros: "<< i-NEventPerRunBruit<<endl;				
				EventNumber = i;
            			Source->GetRandom2(PosXEvent,PosYEvent);
				//indice = 1;            			
				Events->Fill();
            			EventID++;
			
			}

		}
        
    	}
    	fileTree->Write();
        fileTree->Close();
}




double** centrecercle(double R,double r,double XCentreDuCercle,double YCentreDuCercle,double teta0,int N )
{       	
	
	double** x = new double*[N];
	for (int t = 0; t < N ; t++)
	{
		x[t] = new double [2];
	}        
		
	double q = 2*TMath::ATan(r/R);
        for (int i = 0 ; i < N ; i++ )
            {
             x[i][0] = R*TMath::Cos(i*q + teta0 ) + XCentreDuCercle;
             x[i][1] = R*TMath::Sin(i*q + teta0 ) + YCentreDuCercle;
            }	
	return x;
}

double* evenement(double* teta0,double* R,double r,double* RunPositionX , double* RunPositionY,int NrunsValide , int* N,double* XzoneON ,double* YzoneON,double Rpetit,int cb , int g)
{	
	int tour = round(r/Rpetit);        
	int comptNB = 0;
	double** XInterieurCercle = centrecercle(r-Rpetit,Rpetit,XzoneON[0],YzoneON[0],0,round((TMath::Pi()*(r-Rpetit))/Rpetit));
	double norme;	
	double dO;	
	double DZoneRun;
	int NB;
	int c = 0 ;	
	double teta0OFF;	
	int nbOn;	
	int nbOn2;	
	nbOn = round((TMath::Pi()*(r-Rpetit))/0.01);
	nbOn2 = round((TMath::Pi()*(r))/0.01);	
	int poinint = 30;	
	double TabNB[poinint];		
	double** XZoneON = centrecercle(r-Rpetit,0.01,XzoneON[0],YzoneON[0],0,nbOn);	
	double** XZoneON2 = centrecercle(r,0.01,XzoneON[0],YzoneON[0],0,nbOn2);
	for (int i = 0 ; i<poinint ; i++)
	{	
		
		norme = (XInterieurCercle[i][0]) * (XInterieurCercle[i][0])  + (XInterieurCercle[i][1]) * (XInterieurCercle[i][1]); 
	 	dO = 2*TMath::ATan((Rpetit) / (r - Rpetit));		
		DZoneRun = sqrt( (R[0] + ((Rpetit-r) * TMath::Cos(i*dO)) ) * (R[0] + ((Rpetit-r) * TMath::Cos(i*dO)) ) + ((r - Rpetit) * TMath::Sin(i*dO)) * ((r - Rpetit) * TMath::Sin(i*dO)) ); 
		NB = round((TMath::Pi()*DZoneRun)/Rpetit);		
		double** YExterieurCercle = centrecercle(DZoneRun,Rpetit,RunPositionX[0],RunPositionY[0],teta0OFF,NB);						
		comptNB = NB + comptNB;
		for (int p = 0  ; p<NB ; p++) 
		{ 
			if ( (YExterieurCercle[p][0] * YExterieurCercle[p][0]) + (YExterieurCercle[p][1] * YExterieurCercle[p][1]) < r*r ) 
			{
				comptNB = comptNB-1;		
				//cout<< p <<"/"<<NB<<"      "<< i <<endl;			
				NB = NB-1;			
			}	
			else 
			{
		
				//cout<< p <<"/"<<NB<<"      "<< i <<endl;	
				
			}		
		TabNB[i] = NB+1 ; 			
		}
	}	
	double x[comptNB + poinint];
	double y[comptNB + poinint];	
	double* xOn = new double[nbOn];
	double* yOn = new double[nbOn];	
	double* xOn2 = new double[nbOn2];
	double* yOn2 = new double[nbOn2];	
	//cout<< comptNB <<endl;	
	for (int l = 0 ; l<poinint ; l++)
	{	
		norme = (XInterieurCercle[l][0]) * (XInterieurCercle[l][0])  + (XInterieurCercle[l][1]) * (XInterieurCercle[l][1]);  
	 	dO =2*TMath::ATan((Rpetit) / (r - Rpetit));		
		//DZoneRun= sqrt(TMath::Sin(l*dO) * TMath::Sin(l*dO) * norme + (R[0] - sqrt(norme) * TMath::Cos(l*dO)) * (R[0] - sqrt(norme) * TMath::Cos(l*dO)));  
		DZoneRun = sqrt( (R[0] + ((Rpetit-r) * TMath::Cos(l*dO)) ) * (R[0] + ((Rpetit-r) * TMath::Cos(l*dO))) + ((r - Rpetit) * TMath::Sin(l*dO)) * ((r - Rpetit) * TMath::Sin(l*dO)) );
		NB = round((TMath::Pi()*DZoneRun)/Rpetit);					
		teta0OFF = TMath::Pi() - TMath::ATan((TMath::Sin(l*dO) * sqrt(norme)) /(R[0] - (sqrt(norme) * TMath::Cos(l*dO))));  		
		double** YExterieurCercle = centrecercle(DZoneRun,Rpetit,RunPositionX[0],RunPositionY[0],teta0OFF,NB);
		//cout<<"je sui dans la premiere boucle pour le tour: " <<l<<"     " << DZoneRun << "            " <<  teta0OFF << "     " << sqrt(norme) <<"     " << R[0] <<endl;		
				
		for (int p = 0  ; p<NB ; p++) 
		{ 
			if ( (YExterieurCercle[p][0] * YExterieurCercle[p][0]) + (YExterieurCercle[p][1] * YExterieurCercle[p][1]) >= r*r or (p == 0 ) )
			{				
				x[c] = YExterieurCercle[p][0];
				y[c] = YExterieurCercle[p][1];
				c = c + 1 ;		
				//cout<<x[c-1]<<"    "<<x[c-1]<<"    "<< c <<endl;		
			}		
		}	
	
	}
		
	//cout<<nbOn<<endl;	
	for (int t = 0 ; t<nbOn2 ; t++)
	{
		if ( t < nbOn) 		
		{
			xOn[t] = XZoneON[t][0];
			yOn[t] = XZoneON[t][1];
			xOn2[t] = XZoneON2[t][0];
			yOn2[t] = XZoneON2[t][1];		
		}	
		else
		{
			xOn2[t] = XZoneON2[t][0];
			yOn2[t] = XZoneON2[t][1];
		}	
	}	
	DataGenerator(RunPositionX, RunPositionY,NrunsValide,r);
	TFile* fileTree = new TFile("test.root");
        TTree* Events = (TTree*) fileTree->Get("DSTevents");
	int compteur = 0 ;
	double* moyenneoff = new double[poinint] ;
	int n =0 ;
	int ef= 0 ;	
	for (int i = 0 ; i < poinint ; i ++ )
	{
		for (int p = 0 ; p<TabNB[i] ; p ++)
		{	
			
			TString selection = "((PosXEvent-";
			selection += x[compteur];
			selection += ")^2+(PosYEvent-";
			selection += y[compteur];
			selection += ")^2) < ";
			selection += Rpetit;
			selection += "^2";
			compteur = compteur +1;			
			if ( p == 0 ) 
			{
				compteur = compteur +1;		
				
			}

			else
			{
				Events->Draw(">>mylistof_events",selection);
				TEventList *list = (TEventList*)gDirectory->Get("mylistof_events");						
				n = n + list->GetN();  			
				//cout<< list->GetN()<<"     "<<i <<endl;				
				compteur = compteur + 1;
				ef = ef + 1 ;			
			}
 			
		}
		
		moyenneoff[i] = double(n)/ef;		
		ef = 0;
		n = 0;

	} 		
	if (g == cb-1)
	{		
		TGraph* gr = new TGraph(comptNB,x,y);
		//TGraph* gron = new TGraph(nbOn,xOn,yOn);	
		TGraph* gron2 = new TGraph(nbOn2,xOn2,yOn2);	
		gr -> Draw("A*");
		//gron->SetLineColor(2);	
		//gron -> Draw("SAME");
		gron2->SetLineColor(2);	
		gron2 -> Draw("SAME");
	}
	fileTree -> Close();	
	return moyenneoff;
}

void test11()
{ 
    
    
    
    double Rpetit = 0.1;
    int NrunsValide = 0;
    const int Nruns = 1;
    double XzoneON[Nruns] = {0};
    double YzoneON[Nruns] = {0};
    double RunPositionX[Nruns] = {2.5};
    double RunPositionY[Nruns] = {0};
    double R[Nruns];
    int N[Nruns];
    double teta0[Nruns];  
    double r[1] = {2};
    int cb = 1;
    int pp = 70;
    double* exces = new double[cb];
    
    for (int i = 0 ; i < Nruns ; i++)
   {
    	R[i] = sqrt((XzoneON[i]-RunPositionX[i])*(XzoneON[i]-RunPositionX[i])+ (YzoneON[i] - RunPositionY[i])*(YzoneON[i] - RunPositionY[i]));
       	N[i] = round((TMath::Pi()*R[i])/r[0]);
       	if ((YzoneON[i] - RunPositionY[i]) < 0 && (XzoneON[i] - RunPositionX[i]) < 0 )  
	{		
		teta0[i]= TMath::ATan((YzoneON[i] - RunPositionY[i]) / (XzoneON[i] - RunPositionX[i])) + TMath::Pi();
	}	
    	if ((YzoneON[i] - RunPositionY[i]) < 0 && (XzoneON[i] - RunPositionX[i]) > 0 )  
	{
		teta0[i]= TMath::ATan((YzoneON[i] - RunPositionY[i]) / (XzoneON[i] - RunPositionX[i])) ;
	}
      
   	if ((YzoneON[i] - RunPositionY[i]) > 0 && (XzoneON[i] - RunPositionX[i]) > 0 ) 
	{
		teta0[i]= TMath::ATan((YzoneON[i] - RunPositionY[i]) / (XzoneON[i] - RunPositionX[i])) ;
	} 
   
       	if ((YzoneON[i] - RunPositionY[i]) >= 0 && (XzoneON[i] - RunPositionX[i]) < 0 ) 
	{
		teta0[i]= TMath::ATan((YzoneON[i] - RunPositionY[i]) / (XzoneON[i] - RunPositionX[i])) + TMath::Pi();
	}	
        if ((XzoneON[i] - RunPositionX[i]) == 0 && (YzoneON[i] - RunPositionY[i]) > 0 )
	{
	   	teta0[i] = TMath::Pi()/2;
	}
  	if ((XzoneON[i] - RunPositionX[i]) == 0 && (YzoneON[i] - RunPositionY[i]) < 0 )
	{
		teta0[i] = -TMath::Pi()/2;
        } 
   	if ((RunPositionX[i] * RunPositionX[i]) + (RunPositionY[i] * RunPositionY[i]) > (r[0]*r[0]))
	{		 
  		NrunsValide = NrunsValide+1;	
			
	}
  }   
   double* XzoneONNew = new double[NrunsValide];
   double* YzoneONNew = new double[NrunsValide]  ;
   double* RunPositionXNew = new double[NrunsValide];
   double* RunPositionYNew  = new double[NrunsValide];
   double* RNew = new double[NrunsValide];
   double* teta0New = new double[NrunsValide]; 	
   int* NNew = new int[NrunsValide];
   int NombreEcxlus = 0; 
   double moyenne1[cb];
   double moyenne2[cb];
   double  xmin = 0 ;  
   double  ymin = 0 ; 
   double  xmax = 40;
   double  ymax = 6;
   double a ;
   for(int i = 0 ; i< Nruns ; i++)
   {
	if ((RunPositionX[i] * RunPositionX[i]) + (RunPositionY[i] * RunPositionY[i]) > (r[0]*r[0]))
	{	
		XzoneONNew[i-NombreEcxlus] = XzoneON[i]; 
		YzoneONNew[i-NombreEcxlus] = YzoneON[i];
  		RunPositionXNew[i-NombreEcxlus] = RunPositionX[i];
		RunPositionYNew[i-NombreEcxlus] = RunPositionY[i];
		RNew[i-NombreEcxlus] = R[i];
		teta0New[i-NombreEcxlus] = teta0[i];
		NNew[i-NombreEcxlus] = N[i];
 	} 
        else
	{
		cout<<"Le run qui a pour x: "<<RunPositionX[i]<<" et y: "<<RunPositionY[i]<<" a ete ecxlus"<<endl;
		NombreEcxlus = NombreEcxlus + 1 ;
	}
  }  
  for (int g = 0 ; g < cb ; g++)
  { 
	double* moyenne =evenement(teta0New,RNew,r[0],RunPositionXNew,RunPositionYNew,NrunsValide,NNew,XzoneONNew,YzoneONNew,Rpetit , cb, g  );
	moyenne1[g] = moyenne[0];
	moyenne2[g] = moyenne[1];
  	cout<<g<<endl;
  }
  TCanvas* canvascercle = new TCanvas("cercle","cercle");
  TGraph* gr = new TGraph(cb,moyenne1,moyenne2);
  //gr->GetXaxis()->SetRangeUser(0,40);
  //gr->GetYaxis()->SetRangeUser(0,7);
  gr -> Draw("A*");
  a = gr -> GetCorrelationFactor();
  cout<<a<<endl;
}   

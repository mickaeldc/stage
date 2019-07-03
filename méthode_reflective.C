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




void DataGenerator(double* RunPositionX,double* RunPositionY,int Nruns,int FrequenceSource,double r)
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
    	TFile* fileTree = new TFile("fichier_r=1.root","RECREATE");
    	TTree* Events = new TTree("DSTevents","DSTevents");
    
    	Events->Branch("PosXEvent",&PosXEvent,"PosXEvent/D");
    	Events->Branch("PosYEvent",&PosYEvent,"PosYEvent/D");
    	Events->Branch("EventNumber",&EventNumber,"EventNumber/D");
    	Events->Branch("RunNumber",&RunNumber,"RunNumber/D");
	//Events->Branch("indice",&indice,"indice/D");
	
	int FrequenceBruit = 265;    	
	int TimePerRun =  28*60;	
	double NEventPerRunBruit = (FrequenceBruit * TimePerRun)/100;
    	double NEventSource = (FrequenceSource * TimePerRun)/100;
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

double evenement(double* teta0,double* R,double r,double* RunPositionX , double* RunPositionY,int NrunsValide , int* N,int v,int cb , int pp,int w,double* XzoneON ,double* YzoneON,int FrequenceSource)
{	
		
	int c = 0;
	double a = 4;		
	double m = 0.; 	
	int g = 0;	
	double XExclus1 = 0 ;
	double XExclus2 = 0 ;
	double yExclus1 = 0 ;
	double yExclus2 = 0 ;	
	double* TabRecEcX = new double[5];
	double* TabRecEcY = new double[5];	
	int NrunsValideNew =  0;
		
	TabRecEcX[0] = XExclus1;
	TabRecEcX[1] = XExclus2;
	TabRecEcX[2] = XExclus2;
	TabRecEcX[3] = XExclus1;
	TabRecEcX[4] = XExclus1;	
	TabRecEcY[0] = yExclus2;
	TabRecEcY[1] = yExclus2;
	TabRecEcY[2] = yExclus1;
	TabRecEcY[3] = yExclus1;
	TabRecEcY[4] = yExclus2;
        
	int mq = 0;
	int NombreDePoint;
	int CompteurDePoint;
	int* OFFValide = new int[NrunsValide];
	int CompteurDeZone = 0 ; 	
	
	
	for (int p = 0 ; p<NrunsValide ; ++p)
        {		
		
		CompteurDeZone = 0;		
		double** X = centrecercle(R[p],r,RunPositionX[p],RunPositionY[p],teta0[p],N[p]);
		for (int wx = 0 ; wx < N[p] ; wx++)
		{
			CompteurDePoint = 0;			
			NombreDePoint = round((TMath::Pi()*r)/0.01);			
			double** Y = centrecercle(r,0.01,X[wx][0],X[wx][1],0,NombreDePoint);
			for (int gf = 0 ; gf<NombreDePoint ; gf++)
			{
				if ((XzoneON[p] -Y[gf][0] )*(XzoneON[p] -Y[gf][0] ) + (YzoneON[p] - Y[gf][1])*(YzoneON[p] - Y[gf][1]) <= (a*(r*r)))
				{
					CompteurDePoint = CompteurDePoint + 1 ;
				}
				
			}

			if ((CompteurDePoint == 0) or ((X[wx][0] == XzoneON[p]) && (X[wx][p] == YzoneON[p])))
			{
				CompteurDeZone= CompteurDeZone +1;
				g = g + 1;
			}			
		}		
		OFFValide[p] = CompteurDeZone;
		if ((OFFValide[p] != 0) && (OFFValide[p] != 1)) 
		{
			NrunsValideNew = NrunsValideNew+1;	
		}	
	}
	int hg;	
	int* OFFValideNew = new int[NrunsValideNew];	
	double* RNew = new double[NrunsValideNew];
	double* teta0New = new double[NrunsValideNew];
	int* NNew = new int[NrunsValideNew];
	double* XRunsNew = new double[NrunsValideNew];		
    	double* YRunsNew = new double[NrunsValideNew];  			
	
	for (int i = 0 ; i<NrunsValide ; i++)
	{
				
		if ((OFFValide[i] != 0) && (OFFValide[i] != 1)) 
		{		
			OFFValideNew[i-mq] = OFFValide[i] ;			
			XRunsNew[i-mq]= RunPositionX[i];
			YRunsNew[i-mq]= RunPositionY[i];
			teta0New[i-mq] = teta0[i];		
			NNew[i-mq] = N[i];	
			RNew[i-mq] = R[i];		
			//cout<<OFFValideNew[i-mq]<<endl;		
		}
		else 
		{
			mq = mq + 1 ;
			//cout<<"Le run qui a pour x: "<<RunPositionX[i]<<" et y: "<<RunPositionY[i]<<" a ete ecxlus"<<endl;
		}	
	}	
	double* EventON = new double[NrunsValideNew];	
	double* MoyenneOFF = new double[NrunsValideNew];	
	int n;        
	double* exces = new double[NrunsValideNew];	
	double* x =   new double[g];
	double* y =   new double[g];
	int compteur = 0;	
	int ugc =  0 ;		
	DataGenerator(XRunsNew, YRunsNew,NrunsValideNew,FrequenceSource,r);	
	TFile* fileTree = new TFile("fichier_r=1.root");
        TTree* Events = (TTree*) fileTree->Get("DSTevents");  	
	for (int i = 0 ; i < NrunsValideNew ; i++)
	{	
		hg = OFFValideNew[i] ;		
		double** X = centrecercle(RNew[i],r,XRunsNew[i],YRunsNew[i],teta0New[i],NNew[i]);		
		double* XNew = new double[hg];
		double* YNew = new double[hg];
		ugc =  0 ;		   	    		
		for (int wx = 0 ; wx < NNew[i] ; wx++)
		{
			CompteurDePoint = 0;			
			NombreDePoint = round((TMath::Pi()*r)/0.01);			
			double** Y = centrecercle(r,0.01,X[wx][0],X[wx][1],0,NombreDePoint);
			for (int gf = 0 ; gf<NombreDePoint ; gf++)
			{
				if ((XzoneON[i] -Y[gf][0] )*(XzoneON[i] -Y[gf][0] ) + (YzoneON[i] - Y[gf][1])*(YzoneON[i] - Y[gf][1]) <= (a*(r*r)))
				{
					CompteurDePoint = CompteurDePoint + 1 ;
				}
				
			}

			if ((CompteurDePoint == 0) or ((X[wx][0] == XzoneON[i]) && (X[wx][1] == YzoneON[i])))

			{
				XNew[wx-ugc] =X[wx][0];
				YNew[wx-ugc] =X[wx][1];
				//cout<<"ajouter pour le run: "<< i << " pour x: " <<X[wx][0] <<" et pour y: " <<X[wx][1] << " place: "<< (wx-ugc) + 1 <<"/"<<hg<<endl;
			}
			else
			{
				ugc = ugc + 1;
				//cout<<"la zone off qui a pour x: "<<X[wx][0]<<" et pour y: "<< X[wx][1] << " a ete ecxlus pour le run: "<< i<<endl;
			}		
		}			
		NNew[i] = NNew[i] - ugc;		
		if (i !=0)	
		{    
                    	m = m/(compteur);	
			MoyenneOFF[i-1] = m;			
				m = 0;
		}		
		compteur = 0; 	
		for (int t = 0 ; t < hg ; t++)
		{		
			TString selection = "((PosXEvent-";
			selection += XNew[t];
			selection += ")^2+(PosYEvent-";
			selection += YNew[t];
			selection += ")^2) < ";
			selection += r;
			selection += "^2";
			selection += "&& RunNumber == ";
			selection += i;	
		     	if (t==0)
			{ 				
				
				x[c] = XNew[t];
			      	y[c] = YNew[t];				
				c = c + 1;				
				Events->Draw(">>mylistof_events",selection);
				TEventList *list = (TEventList*)gDirectory->Get("mylistof_events");
				n = list->GetN();				
				EventON[i] = n;	
				//cout<< "Il y a "<<n<<" évenements dans la zone ON du "<< i+1 << " ème run!!"<<endl;			
			}			
			
			if((t != 0) && (((X[t][0]<XExclus1) or (X[t][0]>XExclus2)) or ((X[t][1]<yExclus1) or (X[t][1]>yExclus2))))	
			{							
				x[c] = XNew[t];
			        y[c] = YNew[t];				
				compteur = compteur +1;				
				c = c + 1;				
				Events->Draw(">>mylistof_events",selection);
				TEventList *list = (TEventList*)gDirectory->Get("mylistof_events");
				n = list->GetN();				
				m = m + n;						 	
				//cout<< "Il y a "<<n<<" évenements dans la zone OFF numéro "<< t+1 <<"/"<<N[i] << " du "<< i+1<< " ème run!!"<<"pour R = "<<R[i]<<endl;
						
			}		
		} 	
	}      
	
	MoyenneOFF[NrunsValideNew-1] = m/compteur; 	
	double SommeExces = 0 ; 
	double DeltaEcxes = 0;
	double SommeOn = 0;
	double SommeOFF = 0;
		
	for (int u = 0 ; u < NrunsValideNew ; u++)
	{
	    exces[u] = EventON[u] -MoyenneOFF[u];  
	    SommeExces = SommeExces + exces[u]; 
	    SommeOn = SommeOn + EventON[u];
	    SommeOFF = SommeOFF + MoyenneOFF[u];
	    TString nom = "event du run numero";
	    nom += u;
	    //cout<< "la moyenne des évenements dans la zone ONN du run"<< u+1 <<" = " << EventON[u]<<endl;			    
	    //cout<< "la moyenne des évenements dans les zones OFF du run"<< u+1 <<" = " << MoyenneOFF[u]<<endl;
	    //cout<< "la soustraction des évenements dans les zones OFF et on du run"<< u+1 <<" = " << exces[u]<<endl;
	    TString condition  = " RunNumber == ";
            condition += u; 
	}		
	DeltaEcxes = sqrt(SommeOn + SommeOFF);
	delete[] EventON;	
	delete[] MoyenneOFF;	
	cout<<"SommeExces/DeltaEcxes = "<<SommeExces/DeltaEcxes<<endl;
	int nombre = round((TMath::Pi()*sqrt(a)*r)/0.001);	
	double** Y = centrecercle(sqrt(a)*r,0.001,XzoneON[0],YzoneON[0],0,nombre);	
	double* Xexclus = new double[nombre]; ; 
	double* Yexclus = new double[nombre]; 	
	for (int f = 0 ; f < nombre ; f++)
	{
		Xexclus[f] =Y[f][0];  
		Yexclus[f] =Y[f][1];
	}
		
	
	TCanvas* canvascercle = new TCanvas("cercle","cercle");
	TGraph* gr = new TGraph(g,x,y);       
 	gr->Draw("A*");		
	TGraph* grExclus = new TGraph(nombre,Xexclus,Yexclus);	
	grExclus->SetLineColor(2);
	grExclus->Draw("SAME");	
	delete[] Xexclus;
	delete[] Yexclus;	
	delete[] x ;
	delete[] y ;	
	TGraph* grRecExclus = new TGraph(5,TabRecEcX,TabRecEcY);	
	grRecExclus->SetLineColor(2);
	grRecExclus->Draw("SAME");			
	fileTree->Close();	
	delete [] exces;
	return SommeExces/DeltaEcxes;
	
	
}

void test10()
{ 
    
    TFile* file = new TFile("fichiertest.root","RECREATE");
    TTree* cara = new TTree("Caracteristique","Caracteristique");
    double Mean;
    double MeanError;
    double Std;
    double StdError;
    double FrequenceSource = 0;
    
    cara->Branch("Mean",&Mean,"Mean/D");  
    cara->Branch("MeanError",&MeanError,"MeanError/D");
    cara->Branch("Std",&Std,"Std/D");
    cara->Branch("StdError",&StdError,"StdError/D");
    cara->Branch("FrequenceSource",&FrequenceSource,"FrequenceSource/D");
    
    int NrunsValide = 0;
    const int Nruns = 15;
    double XzoneON[Nruns] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double YzoneON[Nruns] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    double RunPositionX[Nruns] = {-3,-1.8,-0.6,0.6,1.8,3.0,-1.8,-0.6,0.6,1.8,-1.8,-1.8,-0.6,0.6,1.8};
    double RunPositionY[Nruns] = {0.8,0.8,0.8,0.8,0.8,0.8,2,2,2,2,2,3.2,3.2,3.2,3.2};
    double R[Nruns];
    int N[Nruns];
    double teta0[Nruns];  
    double r[1] = {0.1};
    int cb = 1;
    int pp = 1;
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
   
   for (int t = 0 ; t< pp ; t++)
   {    		
	TString compteur = "simulation pour f=";
	compteur += FrequenceSource + 1;
	compteur += "hz et r=";
	compteur += r[0];
	compteur +="°";
	cout<<compteur<<endl;		
	FrequenceSource =  FrequenceSource + 1;	
	TH1F* h1 = new TH1F("histo", "histo", 75, -10, 10);        
	for (int h = 0 ; h<cb ; h++)
        {	  	 		
		cout<<"tour numero "<< h <<endl;
	 	
   	 	exces[h] =  evenement(teta0New,RNew,r[0],RunPositionXNew,RunPositionYNew,NrunsValide,NNew,h,cb,pp,t,XzoneONNew,YzoneONNew, FrequenceSource);
         	h1->Fill(exces[h]);
    	} 	
	
	Mean = h1 -> GetMean();
	MeanError = h1 -> GetMeanError();
	Std = h1 -> GetStdDev();
	StdError = h1 -> GetStdDevError();
        cara->Fill();
        h1 -> Delete(); 
   }

  file->Write();
  file->Close();
  
}   

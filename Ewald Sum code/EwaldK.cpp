/*
Writen by Zhao Yihao at 2022/9/1

zhaoyihao@protonmail.com
*/

#include<iostream>
#include<fstream>
#include<cmath>
#include<vector>
#include<string>
#include<sstream>
#include<algorithm>

void ReadCoor( int& nT, std::vector<std::string>& ele, std::vector<std::vector<double>>& r );

void GenerateCoeffi( double& r4pie );

void ForceK( int& nT, double& r4pie, std::vector<std::vector<double>>& r,
              std::vector<std::vector<double>>& L, std::vector<double>& q,
              std::vector<std::vector<double>>& F );

int main()
{

  int nO = 64, nc = 192, nX = 256, nT = 448;
  double lx = 12.4171910, ly = 12.4171910, lz = 12.4171910;

  // Box size
  std::vector<std::vector<double>> L;
  {
    std::vector<double> L1;
    L1.push_back(lx); L1.push_back(1.0/lx);
    L.push_back(L1);
    std::vector<double> L2;
    L2.push_back(ly); L2.push_back(1.0/ly);
    L.push_back(L2);
    std::vector<double> L3;
    L3.push_back(lz); L3.push_back(1.0/lz);
    L.push_back(L3);
  }

  // Atom coor
  std::vector<std::string> ele;
  std::vector<std::vector<double>> r;
  ReadCoor( nT, ele, r );

  // Atom charge
  std::vector<double> q;

  for ( int i = 0; i < nT; i++ )
  {

    if ( ele[i] == "O" ) { q.push_back(6.0); }
    else if ( ele[i] == "H" ) { q.push_back(1.0); }
    else if ( ele[i] == "AX" ) { q.push_back(-2.0);  }
    else if ( ele[i] == "BX" ) { q.push_back(-2.0); }
    else
    {
      std::cout << "there is " << ele.size() <<" ele, not" << nT << std::endl;
    }
     
  }
  
  double r4pie;
  GenerateCoeffi( r4pie );

  std::vector<std::vector<double>> F( nT, std::vector<double>(3, 0) );

  ForceK( nT, r4pie, r, L, q, F );

  std::ofstream out("ForceK_Cpp.xyz");
  out.precision(16);

  for ( int i = 0; i < nT; i++ )
  {
    out << F[i][0] <<"  "<<F[i][1]<<"  "<< F[i][2] << std::endl;
  }
  
  return 0;

}

void ReadCoor( int& nT, std::vector<std::string>& ele, std::vector<std::vector<double>>& r )
{

  std::ifstream coor_in("coor.xyz", std::ios:: in);

  if ( !coor_in.is_open() )
  {
    std::cout << "File is not exist" << std::endl;
  }

  double x, y, z;
  std::string hang;
  std::string ele_lin;

  getline(coor_in, hang);
  getline(coor_in, hang);

  for ( int i = 0; i < nT; i++ )
  {

    getline(coor_in, hang);

    if ( coor_in.eof() )
    {
     std::cout << "end of file" << std::endl;
    }

    std::istringstream is(hang);

    is >> ele_lin >> x >> y >> z;

    ele.push_back(ele_lin);

    std::vector<double> r_lin;
    r_lin.push_back(x);
    r_lin.push_back(y);
    r_lin.push_back(z);

    r.push_back(r_lin);

  }

return;

}

void GenerateCoeffi( double& r4pie )
{
  double pi = acos(-1.0);
  double ele_to_clbp19 = 1.6021766340;
  double epslonp12 = 8.85418781280;
  double nan23=6.022140760;

  r4pie = pow(ele_to_clbp19,2)*nan23/(4.0*pi*epslonp12)*1.0e6; // 10J/mol/ans

return;

}


void ForceK( int& nT, double& r4pie, std::vector<std::vector<double>>& r,
              std::vector<std::vector<double>>& L, std::vector<double>& q, 
              std::vector<std::vector<double>>& F )
{

/*

 nT is the total number of particles;
 r4pi4 is the Unit coefficient of conversion;
 r is an nT*3 matrix storing the coordinates of the particles;
 L[][0] is ( lx, ly, lz ), l[][1] is (1/lx, 1/ly, 1/lz );
 q is a vector of size nT that stores the charge of the particle;
 F is an nT*3 matrix that returns the force on the particle; 
 
*/

  double alpha = 1.0/4.50, pi = acos(-1.0);
  double alpha2 = alpha*alpha, twopi = 2.0*pi;
  double p4alp2 = -0.250/alpha2;
  int kLmax = 10, kMmax = 10, kNmax = 10;
  double rcpcut = std::min( {double(kLmax)*L[0][1], double(kMmax)*L[1][1], double(kNmax)*L[2][1]} );
  rcpcut = rcpcut*1.050*twopi;
  double rcpct2 = rcpcut*rcpcut;

  std::vector<std::vector<double>> eLc(nT, std::vector<double>(kLmax+1) ), 
                                   eLs(nT, std::vector<double>(kLmax+1) );
  std::vector<std::vector<double>> eMc(nT, std::vector<double>(kMmax+1) ),
                                   eMs(nT, std::vector<double>(kMmax+1) );
  std::vector<std::vector<double>> eNc(nT, std::vector<double>(kNmax+1) ), 
                                   eNs(nT, std::vector<double>(kNmax+1) );

  for ( int i = 0; i < nT; i++ )
  {

    eLc[i][0] = 1.0;  eLs[i][0] = 0;
    eMc[i][0] = 1.0;  eMs[i][0] = 0;
    eNc[i][0] = 1.0;  eNs[i][0] = 0;

    double ssX = r[i][0]*L[0][1];
    double ssY = r[i][1]*L[1][1];
    double ssZ = r[i][2]*L[2][1];

    eLc[i][1] = cos(twopi*ssX); eLs[i][1] = sin(twopi*ssX);
    eMc[i][1] = cos(twopi*ssY); eMs[i][1] = sin(twopi*ssY);
    eNc[i][1] = cos(twopi*ssZ); eNs[i][1] = sin(twopi*ssZ);

  } 


  for ( int kM = 2; kM <= kMmax; kM++ )  
  {
    for ( int j = 0; j < nT; j++ )
    {
      eMc[j][kM] = eMc[j][kM-1]*eMc[j][1] - eMs[j][kM-1]*eMs[j][1];
      eMs[j][kM] = eMs[j][kM-1]*eMc[j][1] + eMc[j][kM-1]*eMs[j][1];
    }
  }

  for ( int kN = 2; kN <= kNmax; kN++ )
  {
    for ( int j = 0; j < nT; j++ )
    {
      eNc[j][kN] = eNc[j][kN-1]*eNc[j][1] - eNs[j][kN-1]*eNs[j][1];
      eNs[j][kN] = eNs[j][kN-1]*eNc[j][1] + eNc[j][kN-1]*eNs[j][1];
    }
  }
     
////////////////////////////////////////////////////////////////////


  int kMmin = 0, kNmin = 1;

  for ( int LL = 0; LL <= kLmax; LL++ )
  {
    int kL = LL;
    double rkx = twopi*double(kL)*L[0][1];

    if ( kL == 1 ) 
    {
      for ( int j = 0; j < nT; j++ )
      {
        eLc[j][0] = eLc[j][1]; 
        eLs[j][0] = eLs[j][1];
      }
    }
    else if ( kL > 1 ) 
    {
      for ( int j = 0; j < nT; j++ )
      {
        double cs = eLc[j][0];
        eLc[j][0] = cs       *eLc[j][1] - eLs[j][0]*eLs[j][1];
        eLs[j][0] = eLs[j][0]*eLc[j][1] + cs       *eLs[j][1];
      }
    }
 
  
    for ( int MM = kMmin; MM <= kMmax; MM++ )
    {

      int kM = abs(MM);
      double rky = twopi*double(MM)*L[1][1];
      std::vector<double> clm(nT), slm(nT);
      
      if ( MM >= 0 ) 
      {
        for ( int j = 0; j < nT; j++ )
        {
          clm[j] = eLc[j][0]*eMc[j][kM] - eLs[j][0]*eMs[j][kM];
          slm[j] = eLs[j][0]*eMc[j][kM] + eLc[j][0]*eMs[j][kM];
        }
        
      }
      else 
      {
        for ( int j = 0; j < nT; j++ )
        {
          clm[j] = eLc[j][0]*eMc[j][kM] + eLs[j][0]*eMs[j][kM];
          slm[j] = eLs[j][0]*eMc[j][kM] - eLc[j][0]*eMs[j][kM]; 
        }

      }
 
      for ( int NN = kNmin; NN <= kNmax; NN++ )
      {
        int kN = abs(NN);
        double rkz = twopi*double(NN)*L[2][1];
        double rksq = rkx*rkx + rky*rky + rkz*rkz;
        std::vector<double> coskr(nT), sinkr(nT);
        
        if ( rksq <= rcpct2 ) 
        {
          if ( NN >= 0 ) 
          {
            for ( int j = 0; j < nT; j++ )
            {
              coskr[j] = clm[j]*eNc[j][kN] - slm[j]*eNs[j][kN];
              sinkr[j] = slm[j]*eNc[j][kN] + clm[j]*eNs[j][kN];
            }
          }
          else 
          {
            for ( int j = 0; j < nT; j++ )
            {
              coskr[j] = clm[j]*eNc[j][kN] + slm[j]*eNs[j][kN];
              sinkr[j] = slm[j]*eNc[j][kN] - clm[j]*eNs[j][kN];
            }

          }
         
          double sumqcoskr = 0;
          double sumqsinkr = 0;
  
          for ( int j = 0; j < nT; j++ )
          {
            sumqcoskr = sumqcoskr + q[j]*coskr[j];
            sumqsinkr = sumqsinkr + q[j]*sinkr[j];
          }
     
          double akk = exp(p4alp2*rksq)/rksq;     
          
          for ( int j = 0; j < nT; j++ )
          {
            double FF = akk*q[j]*( sinkr[j]*sumqcoskr - coskr[j]*sumqsinkr );
             
            F[j][0] = F[j][0] + rkx*FF;
            F[j][1] = F[j][1] + rky*FF;
            F[j][2] = F[j][2] + rkz*FF;
          }     
          
        } // if rksq <= rcpct2

      } // int NN 
 
      kNmin = -kNmax;

    } // int MM

    kMmin = -kMmax;

  } // int LL
 

  double Fcorffi = 4.0*twopi*r4pie*L[0][1]*L[1][1]*L[2][1];

  for ( int i = 0; i <= 2; i++ )
  {
    for ( int j = 0; j < nT; j++ )
    {
      F[j][i] = F[j][i]*Fcorffi;
    }
  }

return;

}




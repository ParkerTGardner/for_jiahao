const int   trackbin                    =5 ;
const int   ktbin = 8;
const int ptbin = 6;
const float ptbinbounds_lo[ptbin] = {0.3,0.7,1.0,1.5,2.0,2.5};
const float ptbinbounds_hi[ptbin] = {0.7,1.0,1.5,2.0,2.5,3.0};
const int mm = 1000;
const int   trackbinbounds[trackbin]         = {0, 0, 30, 0, 60};
const int   trackbinboundsUpper[trackbin]    = {mm,30,60, 60,mm};
const float ktbinbounds_lo[ktbin] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
const float ktbinbounds_hi[ktbin] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};
const double        PI = 3.14159265359;
const float         EtaBW = 0.3;
const float         PhiBW = TMath::Pi()/16;
const int bin_WRTJ_Eta        = 150;
const int low_WRTJ_Eta_Bin  =  0;
const int high_WRTJ_Eta_Bin = 10;
const int bin_Eta             = 60;
const int low_Eta_Bin       = -3;
const int high_Eta_Bin      = 3;

const int EPD_xb  = 300;
const int EPD_yb  = 120;
const int EPD_xhi = 10;
const int EPD_xlo = -10;
const int EPD_yhi = 4;
const int EPD_ylo = -4;

const double bin2pi =2*TMath::Pi();
const int bin0 =0;
const int bin1 =1;
const int bin2 =2;
const int bin3 =3;
const int bin4 =4;
const int bin10 =10;
const int bin20 =20;
const int bin30 =30;
const int bin60 =60;
const int bin100 =100;
const int bin200 =200;
const int bin120 =150;
const int bin360 =450;
const int bin150 =150;
const int bin300 =300;
const int bin3000 =3000;
const int bin8000 =8000;
const int bin270=270;

const int bin500=500;
const int bin50=50;
const int bin22=22;


//---------------------------------------------------------------------CUTS
//const float EtaCut      = 0.0;
const float jetEtaCut   = 1.6;
const float jetPtCut_Event    = 550.0;
const float jetPtCut_Jet    = 550.0;

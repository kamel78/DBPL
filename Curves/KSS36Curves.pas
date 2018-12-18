unit KSS36Curves;

interface

uses  Vcl.ComCtrls, System.SysUtils,VLargeIntegers,LargeIntegers,Fp16Arithmetic,Fp8Arithmetic,Fp6Arithmetic,
      Fp2Arithmetic,GeneralTypes,System.classes,VCL.dialogs,ECCFp4,ECCFp6,ECCFp,Fp36Arithmetic,TreeView1;


Type
   GTKSS36=Fp36Int;
   G2KSS36=Fp6Point;
   G1KSS36=FpPoint;
   LargeInt=Lint;
   KSS36CurvesPairingAlgos=(KSS36OptAte);
   StandardKSS36Curves=(scKSS36at256_1,scKSS36at256_2,scKSS36at192_1,scKSS36at192_2);

   KSS36CurvesParamsDefinition=record
                            SecurityLevel:String;
                            u:String;  // the paramater of generation for the KSS36 curve
                            Beta:Integer; // non-square elements of the irredictible polynomial on FP2
                            Sigma:string; // non-square non-cube elements of the irredictible polynomial on FP4
                            Gamma:string; // non-square non-cube elements of the irredictible polynomial on FP8
                            A,B:String; // parametres of the curve   y^2=x^3+Ax+B
                            TwistMode:TTwistModel;
                            GenratorX:String;
                            TwistGeneratorSeed:Word;
                            end;
   ListOfKSS36Params=array of KSS36CurvesParamsDefinition;

 TKSS36Curve=class (TCurve)
                      {Definition of a Curve with respect the the Weistrass model   Y^2=X^3+A*X+B (mod P)}
              private
              _A,_B:Lint;
              _u:Lint;
              _LoopMode:LoopPoweringMode;                     //  Loopmode :Negative/Binary Representation
              _Fp36PoweringMode:FpPoweringMode;             // Final Exponentiation Powering Mode : Beuchat Approach, Karabina Approach, Classical Square and Multiply
              _PairingAlgorithm:KSS36CurvesPairingAlgos;         //  Pairing Algorithme OptAte,R-Ate.....
              power:Lint;
              decomp:LIntArray;
              params:StandardKSS36Curves;
              procedure SetA(Value:String);
              function getA:String;
              procedure SetB(Value:String);
              function getB:String;
              procedure Setu(Value:String);
              function getu:String;
              function getRtw:String;
              function getAtw:String;
              function getBtw:String;
              procedure SetSigma(Value:String);
              function getSigma:String;
              function getGamma:String;
              procedure SetGamma(Value:String);
              procedure SetBeta(Value:integer);
              function getBeta:integer;
              function getR:string;
              function getN:string;
              function getTr:string;
              function getLp:string;
              function getHtw:string;
              function getH:string;
              function getP:string;
              procedure setCoord(value:CoordinatesSystem);
              function getCoord:CoordinatesSystem;
              function GetTwm:TTwistModel;
              procedure RecomputeParametres;
              Procedure FinalPowerOptimalAteKSS1(f:Fp36Int;u:LInt;PoweringMode:FpPoweringMode;var Result:Fp36Int); // Compute f^((p^24-1)/r)
              Function OptAtePairing(Pt:FpPoint;Qt:Fp6Point):Fp36Int;  // Compute Optimal Ate Pairings
              procedure SetLoopMode(Value:LoopPoweringMode);
              procedure SetCurve(Value:TKSS36Curve);
              procedure SetPArams(const Value: StandardKSS36Curves);
            public
              SecurityLevel:String;    /// A symbolic identifier of the curve
              PerferedPoweringMode:LoopPoweringMode;
              Loop:PLIntArrayForm;
              property CoordinatesSystem:CoordinatesSystem read getCoord write SetCoord;
              property Beta:integer read getBeta write SetBeta;
              property Sigma:string read Getsigma write SetSigma;
              property Gamma:string read GetGamma Write SetGamma;
              property u :string read Getu write Setu;
              property A:string read getA write setA;
              property B:string read getB write setB;
              property Atw :string read GetAtw;
              property Btw :string read GetBtw;
              property P:string read GetP;
              property H:string read GetH;
              property N:string read getN;
              property R:String read getR;
              property Rtw:string read getRtw;
              property Htw:string read getHtw;
              property Tr:string read getTr;
              property Lp:string read getLp;
              property TwistMode:TTwistModel read GetTwm;
              Constructor Create(AOwner : TComponent); override;
              destructor Destroy;
              procedure SetStandardCurveParametres(inParametres:StandardKSS36Curves);
              procedure SetCustomCurveParametres(inParametres:KSS36CurvesParamsDefinition);
              procedure GenerateParamsTreeView(Tree:TTreeView);override;
              function Paire(P:FpPoint;Q:Fp6Point):Fp36Int;
              function GetRandomG1Point:FpPoint;
              function GetRandomG2Point:Fp6Point;
              function GetDefautG1Generator:FpPoint;
              function GetDefautG2Generator:Fp6Point;
              function HashToG1Point(id:string):FpPoint;
              function HashToG2Point(id:String):Fp6Point;
              published
              property Fp36PoweringMode:FpPoweringMode read _Fp36PoweringMode write _Fp36PoweringMode;
              property PairingAlgorithm:KSS36CurvesPairingAlgos read _PairingAlgorithm write _PairingAlgorithm ;
              property LoopMode:LoopPoweringMode read _LoopMode write SetLoopMode;
              property Parametres:StandardKSS36Curves read params write Setparams;

            end;

const

 KSS36_256_1_Params:KSS36CurvesParamsDefinition=(SecurityLevel:'256bit';u:'0x61200013330A';Beta:-2;sigma:'0+u*1';A:'0';B:'2';TwistMode:twDType;GenratorX:'-1';TwistGeneratorSeed:55);
 KSS36_256_2_Params:KSS36CurvesParamsDefinition=(SecurityLevel:'256bit';u:'0x126CA209DE61A';Beta:-2;sigma:'0+u*1';A:'0';B:'2';TwistMode:twDType;GenratorX:'-1';TwistGeneratorSeed:55);
 KSS36_192_1_Params:KSS36CurvesParamsDefinition=(SecurityLevel:'192bit';u:'0xC2323FC6A';Beta:-2;sigma:'0+u*1';A:'0';B:'2';TwistMode:twDType;GenratorX:'-1';TwistGeneratorSeed:55);
 KSS36_192_2_Params:KSS36CurvesParamsDefinition=(SecurityLevel:'192bit';u:'0xc2731aa5';Beta:-2;sigma:'0+u*1';A:'0';B:'8';TwistMode:twDType;GenratorX:'1';TwistGeneratorSeed:55);
 KSS36ParamsList:array[0..3] of string=('KSS36_192_1_Params','KSS36_192_2_Params','KSS36_256_1_Params','KSS36_256_2_Params');
 KSS36ImplementedPairingAlgos:array[0..0] of string=('Optimal Ate Pairing');

  procedure ComputeKSS36Parametres(Params: KSS36CurvesParamsDefinition; var Result:PtrCurveParams);
  procedure Register;

implementation

var t:array[0..30] of FP36Int;



procedure ComputeKSS36Parametres(Params: KSS36CurvesParamsDefinition; var Result:PtrCurveParams);
var tmp,tmp1:LInt;
    i:integer;
    m2:HLInt;
    e:Lint;
    _P,_F,_tau:HLInt;
    twtmp:Fp6Int;
    twtmp1:Fp6Int;
    TwoinFp2:Fp2int;
    Genrator,Gtmp:Fp6Point;
    KSSTwistBasePointX,KSSTwistBasePointY:Fp6Int;
begin
if Result=nil then begin
                   new(Result);
                   new(Result.LoopBin);
                   new(Result.LoopNaf);
                   new(Result.LoopTateNaf);
                   new(Result.LoopTateBin);
                   end;
with Result^ do begin
                SecurityLevel:=Params.SecurityLevel;
                Family:=cfKSS36;
                if (Pos('0x',Params.u)<>0)or(Pos('$',Params.u)<>0) then u:=_hex_To_LInt(Params.u)
                else u:=_str_To_LInt(Params.u);

                Tr:=(2*u.Pow(7) + 757*u + 259)/259;
                r:=(u.Pow(12) + 683*u.Pow(6) + 117649)/161061481;
                p:=(u.Pow(14)-4*u.pow(13)+7*u.Pow(12)+683*u.Pow(8)-2510*u.Pow(7) +4781*u.Pow(6)+117649*u.Sqr-386569*u + 823543)/28749;
                n:=p+1-tr;
                h:=n/r;
                if ((not IsLIntPrime(p))or
                (not IsLIntPrime(r)))or
                (n mod r <>0) then
                    raise Exception.Create('Paramétre t invalide pour la construction de la courbe KSS...')
                else begin
                     p.Limit:=2*p.Data.i16[-2];
                     new(p.InverseFordivision);
                     p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
                     p.id:='Y';
                     Lp:=(u).Absolute;
                     end;

                A:=_Str_To_LInt(Params.A);
                B:=_Str_To_LInt(Params.B);
                FieldParam:=InitFieldParams(Params.Beta,p);
                //computing number of points on E(FP6)with cofactor
                e:=params.beta;
                if e.IsASqrMod(P,FieldParam.MontgomeryData) then raise Exception.Create(inttostr(Params.beta)+' is a quadratic Residue modulo P......');
                _tau:=tr.Pow(6)-6*p*tr.Pow(4)+9*p.Sqr*tr.Sqr-2*p.Sqr*p;
                _P:=p.Sqr*p.Sqr*p.Sqr;
                _F:=(4*_P-_Tau.Sqr)/3;
                _Sqrt_HLint(_F,_F);
                m2:=_P+1-(3*_F+_tau)/2;
                _VHCopy_LInt(m2,tmp1);

                Htw:=tmp1/r;    // Twist Cofactor
                Rtw:=r;         ///  The true order is m2=r*htw , so we can choose a subgourp of order R (the cofactor is htw) (https://eprint.iacr.org/2005/133.pdf)

                New(TowerParam2);
                New(TowerParam2.Tower8paprms);
                TowerParam2.Tower8paprms^.FieldParam:=FieldParam;
                TowerParam2.Tower8paprms^.Sigma.SetFieldParams(TowerParam2.Tower8paprms^.FieldParam);
                TowerParam2.Tower8paprms.Sigma.SetFormString(Params.Sigma);
                     ///   should test if beta and sigma are valide parametres
                TwoinFp2.SetFieldParams(FieldParam);
                _Pow_FP2(TowerParam2.Tower8paprms.Sigma,(Result.P.Sqr-1)/2,TwoinFp2);
                if (TwoinFp2.a=1)and(TwoinFp2.b=0) then raise Exception.Create('Paramétre t invalide pour la construction de la courbe KSS...Sigma is a square');
                _Pow_FP2(TowerParam2.Tower8paprms.Sigma,(Result.P.Sqr-1)/3,TwoinFp2);
                if (TwoinFp2.a=1)and(TwoinFp2.b=0) then raise Exception.Create('Paramétre t invalide pour la construction de la courbe KSS...Sigma is a Cube');
                    ///
                TowerParam2.Tower8paprms.pmod8:=p mod 8;
                TwoinFp2.SetFieldParams(TowerParam2.Tower8paprms^.FieldParam);
                TwoinFp2.SetFormString('2+u*0');       // for squaring over FP6
                TowerParam2.Tower8paprms.inv_2.SetFieldParams(TowerParam2.Tower8paprms.FieldParam);
                TowerParam2.Tower8paprms.inv_2:=TwoinFP2.Inverse;
                TowerParam:=TowerParam2.Tower8paprms;

                TowerParam2.pmod8:=TowerParam2.Tower8paprms.pmod8;
                TowerParam2.Sigma:=TowerParam2.Tower8paprms.Sigma;
                TowerParam2.FieldParam:=FieldParam;
                TowerParam2.Gamma2.SetTowerParams(TowerParam2.Tower8paprms);
                TowerParam2.Gamma2.a.SetToZero;
                TowerParam2.Gamma2.b.SetToOne;
                TowerParam2.Gamma2.c.SetToZero;

                BtwFp6.SetTowerParams (TowerParam);
                AtwFp6.SetTowerParams(TowerParam);
                AtwFp6.a.SetToZero;
                AtwFp6.b.SetToZero;
                AtwFp6.c.SetToZero;

                TwistMode:=Params.TwistMode;
                if TwistMode=twMType then _Mul_FP_FP6(B,TowerParam2.Gamma2,BtwFp6) else
                _Mul_FP_FP6(B,TowerParam2.Gamma2.Inverse,BtwFp6);

                LoopBin^:=Lp.ToIntArray;
                LoopNaf^:=Lp.ToNafArray;
                LoopTateBin^:=n.ToIntArray;
                LoopTateNaf^:=n.ToNafArray;
                if (Pos('0x',Params.GenratorX)<>0)or(Pos('$',Params.GenratorX)<>0) then BasePointX:=_hex_To_LInt(Params.GenratorX)
                else BasePointX:=_Str_To_LInt(Params.GenratorX);

                ///   Constructing the base point of the generator for G1
                _Sqr_LInt(BasePointX,tmp);
                _Add_LInt(tmp,A,tmp);
                _Mul_LInt(tmp,BasePointX,tmp1);
                _Add_LInt(tmp1,B,tmp);
                _Mod_LInt(tmp,P,tmp);
                if tmp.IsASqrMod(P,FieldParam.MontgomeryData)
                                              then begin
                                                   ModSquareLInt(tmp,P, BasePointY,FieldParam.MontgomeryData,false);
                                                   if P-BasePointY<BasePointY then BasePointY:=P-BasePointY;
                                                   end
                else raise Exception.Create('Invalide parameter for the generator of the curve');
                ///   Constructing the base point for G2 Generator
                KSSTwistBasePointX.SetTowerParams(TowerParam);
                KSSTwistBasePointY.SetTowerParams(TowerParam);
                KSSTwistBasePointX.a.a.SetToRandom(p,params.TwistGeneratorSeed);
                KSSTwistBasePointX.a.b.SetToRandom(p,params.TwistGeneratorSeed);
                KSSTwistBasePointX.b.a.SetToRandom(p,params.TwistGeneratorSeed);
                KSSTwistBasePointX.b.b.SetToRandom(p,params.TwistGeneratorSeed);
                KSSTwistBasePointX.c.a.SetToRandom(p,params.TwistGeneratorSeed);
                KSSTwistBasePointX.c.b.SetToRandom(p,params.TwistGeneratorSeed);


                repeat
                _Inc_LInt(KSSTwistBasePointX.a.a,1);
                _Sqr_FP6(KSSTwistBasePointX,twtmp);
                _Add_FP6(twtmp,AtwFp6,twtmp);
                _Mul_FP6(twtmp,KSSTwistBasePointX,twtmp1);
                _Add_FP6(twtmp1,BtwFp6,twtmp);
                until twtmp.IsASquare;
                _Sqrt_FP6(twtmp,KSSTwistBasePointY);
                          /// constants for FP2 powering
                TowerParam^.FrobeniusP_Const[0].SetFieldParams(TowerParam^.FieldParam);
                TowerParam^.FrobeniusP_Const[0]:=TowerParam^.Sigma.Pow((P-1)/3);
                for i:=1 to 4 do begin
                                 TowerParam^.FrobeniusP_Const[i].SetFieldParams(TowerParam^.FieldParam);
                                 TowerParam^.FrobeniusP_Const[i]:=TowerParam^.FrobeniusP_Const[i-1]*TowerParam^.FrobeniusP_Const[0];
                                 end;
                          /// constants for FP6 powering
                TowerParam2^.FrobeniusP_Const_Fp6[0].SetTowerParams(TowerParam);
                TowerParam2^.FrobeniusP_Const_Fp6[0]:=TowerParam2^.Gamma2.Pow((P-1)/6);
                for i:=1 to 4 do begin
                                 TowerParam2^.FrobeniusP_Const_Fp6[i].SetTowerParams(TowerParam);
                                 TowerParam2^.FrobeniusP_Const_Fp6[i]:=TowerParam2^.FrobeniusP_Const_Fp6[i-1]*TowerParam2^.FrobeniusP_Const_Fp6[0];
                                 end;
                FrobeniusMapConstX_Fp4.SetTowerParams(TowerParam);
                case TwistMode of
                    twDType: begin
                             FrobeniusMapConstX_Fp6:=TowerParam2^.Gamma2.Pow((p-1)/2);
                             FrobeniusMapConstY_Fp6:=TowerParam2^.Gamma2.Pow(3*(p-1)/4);
                             end;
                twMType:begin
                        FrobeniusMapConstX_Fp6:=(TowerParam2^.Gamma2.Pow((p-1)/2)).Inverse;///
                        FrobeniusMapConstY_Fp6:=(TowerParam2^.Gamma2.Pow(3*(p-1)/4)).Inverse;   //////
                        end;
                end;
                TowerParam.FrobeniusPi3xSigmaSqr:=(TowerParam2^.Gamma2.Pow((p-3)).a);
                TowerParam.FrobeniusPi3xSigmaSqr:=TowerParam.FrobeniusPi3xSigmaSqr*TowerParam2.Sigma;
                TowerParam2.FrobeniusPi3xSigmaSqr:=TowerParam.FrobeniusPi3xSigmaSqr;
                Genrator.SetCurveParams(Result);
                Gtmp.SetCurveParams(Result);
                Gtmp.X:=KSSTwistBasePointX;
                Gtmp.Y:=KSSTwistBasePointY;
                Gtmp.Z.a.SetFormString('1');
                Gtmp.Z.b.SetFormString('0');
                Gtmp.Infinity:=false;
                _Mul_Fp6_FpPoint(htw,Gtmp,Genrator);
                KSS36TwistGeneratorX:=Genrator.X;
                KSS36TwistGeneratorY:=Genrator.Y;
                end;
end;


{ TKSS36Curve }

procedure TKSS36Curve.GenerateParamsTreeView(Tree:TTreeView);
var s:string;
    r:real;
    G:FpPoint;Gtw:Fp6Point;
    tmp:TTreeNode;
begin
Tree.Items.Clear;
tmp:=Tree.Items.Add(nil,'Curve');
Tree.Items.AddChild(tmp,'Family : KSS36');
if CurveParams.A<0 then  s:='y^2=x^3-'+CurveParams.B.Absolute.ToDecimalString
else s:='y^2=x^3+'+CurveParams.B.ToDecimalString;
Tree.Items.AddChild(tmp,'Equation : '+s);
tmp:=Tree.Items.Add(nil,'Initial Parameter x0');
Tree.Items.AddChild(tmp,'x0= '+curveparams.u.ToHexString );
Tree.Items.AddChild(tmp,'Hamming Weight (Bin) ='+IntToStr(HammingWeight(CurveParams.u,False)));
Tree.Items.AddChild(tmp,'Hamming Weight (Neg) ='+IntToStr(HammingWeight(CurveParams.u,True)));
tmp:=Tree.Items.Add(nil,'Field Fp');
Tree.Items.AddChild(tmp,'Prime P = '+CurveParams.P.ToHexString );
Tree.Items.AddChild(tmp,'Size of Fp : '+inttostr(CurveParams.P.BitLength)+' bit');
tmp:=Tree.Items.Add(nil,'Field E(Fp)');
Tree.Items.AddChild(tmp,'Order (#E(Fp)) N = '+CurveParams.N.ToHexString );
Tree.Items.AddChild(tmp,'Frobenius Trace Tr : '+CurveParams.Tr.ToHexString);
tmp:=Tree.Items.Add(nil,'Base Fields G1 : E(Fp)[R]');
Tree.Items.AddChild(tmp,'Order of G1 R = '+CurveParams.R.ToHexString);
Tree.Items.AddChild(tmp,'Size of G1 : '+inttostr(CurveParams.R.BitLength)+' bit');
Tree.Items.AddChild(tmp,'Cofactor of G1 H = '+CurveParams.H.ToHexString);
tmp:=Tree.Items.Add(nil,'Twist');
Tree.Items.AddChild(tmp,'Equation of the Twist : '+'y^2=x^3+'+Curveparams.BtwFp6.toHexString);
Tree.Items.AddChild(tmp,'Degree of the Twist : Sextic');
if Curveparams.TwistMode=twMType then S:='M-Type' else S:='D-Type';
Tree.Items.AddChild(tmp,'Type of the Twist : '+S);
tmp:=Tree.Items.Add(nil,'Field G2 (Twist) :E''(Fp6)[R]');
Tree.Items.AddChild(tmp,'Order of G2 R = '+CurveParams.Rtw.ToHexString);
Tree.Items.AddChild(tmp,'Size of G2 : '+inttostr(CurveParams.Rtw.BitLength)+' bit');
Tree.Items.AddChild(tmp,'Cofactor of G2 Htw = '+CurveParams.Htw.ToHexString);
tmp:=Tree.Items.Add(nil,'Targted Field (GT):');
Tree.Items.AddChild(tmp,'Extension Field : Fp36');
Tree.Items.AddChild(tmp,'Size of GT : '+inttostr(36*CurveParams.P.BitLength)+' bit');
tmp:=Tree.Items.Add(nil,'Security Level');
Tree.Items.AddChild(tmp,'Security in G1 : '+inttostr(CurveParams.R.BitLength div 2)+' bit').ImageIndex:=2;
Tree.Items.AddChild(tmp,'Security in G2 : '+inttostr(CurveParams.R.BitLength div 2)+' bit').ImageIndex:=2;
r:=1.92* exp((1/3)*ln(36*CurveParams.P.BitLength))*sqr(exp((1/3)*(ln(ln(36*CurveParams.P.BitLength)/ln(2)))));
Tree.Items.AddChild(tmp,'Security in GT (GNFS) : '+inttostr(Round(r))+' bit').ImageIndex:=2;
tmp:=Tree.Items.Add(nil,'Implemented Pairings');
Tree.Items.AddChild(tmp,'Otp-Ate');
tmp:=Tree.Items.Add(nil,'Tower Construction');
if Curveparams.FieldParam.Beta<0 then Tree.Items.AddChild(tmp,'Fp2<u>=ExstensionField<u,|u^2+'+inttostr(Curveparams.FieldParam.Beta)+'>')
else Tree.Items.AddChild(tmp,'Fp2<u>=ExstensionField<u,|u^2-'+inttostr(Curveparams.FieldParam.Beta)+'>');
Tree.Items.AddChild(tmp,'Fp6<v>=ExstensionField<v,|v^3-('+Curveparams.TowerParam2.Sigma.toHexString+')>');
Tree.Items.AddChild(tmp,'Fp12<z>=ExstensionField<z,|z^2-w>');
Tree.Items.AddChild(tmp,'Fp36<t>=ExstensionField<t,|t^3-z>');
tmp:=Tree.Items.Add(nil,' = '+floattostrf(CurveParams.P.BitLength/CurveParams.R.BitLength,ffGeneral, 2, 4 ));
tmp:=Tree.Items.Add(nil,'Default G1 Generator');
G.SetCurveParams (CurveParams);
G.SetToDefaultGenerator;
Gtw.SetCurveParams (CurveParams);
Gtw.SetToDefaultGenerator;
Tree.Items.AddChild(tmp,G.toHexString);
tmp:=Tree.Items.Add(nil,'Default G2 Generator');
Tree.Items.AddChild(tmp,Gtw.toHexString);
Tree.FullExpand;
Tree.Selected:=Tree.Items[0];
//Tree.SetFocus;
end;


constructor TKSS36Curve.Create;
begin
inherited Create(AOwner);
CurveParams:=nil;
SetStandardCurveParametres(scKSS36at192_1);
CoordinatesSystem:=csProjective;
_LoopMode:=lpmAuto;
PairingAlgorithm:=KSS36OptAte;
end;

{*******************************************************************************}
destructor TKSS36Curve.Destroy;
begin
Dispose(CurveParams);
inherited;
end;

{*******************************************************************************}
function TKSS36Curve.getA: String;
begin
Result:=_A.ToHexString;
end;

{*******************************************************************************}
function TKSS36Curve.getAtw: String;
begin
Result:=CurveParams^.AtwFp6.toHexString;
end;

{*******************************************************************************}
function TKSS36Curve.getB: String;
begin
Result:=_B.ToHexString;
end;

{*******************************************************************************}
function TKSS36Curve.getBeta:integer;
begin
Result:=CurveParams^.FieldParam.Beta;
end;

{*******************************************************************************}
function TKSS36Curve.getBtw: String;
begin
Result:=CurveParams^.BtwFp6.toHexString;
end;

{*******************************************************************************}
function TKSS36Curve.getCoord: CoordinatesSystem;
begin
Result:=CurveParams^.CoordSys;
end;

function TKSS36Curve.GetDefautG1Generator: FpPoint;
begin
Result.SetToDefaultGenerator;
end;

function TKSS36Curve.GetDefautG2Generator: Fp6Point;
begin
Result.SetToDefaultGenerator;
end;

{*******************************************************************************}
function TKSS36Curve.getGamma: String;
begin
if  TwistMode=twDType then Result:=CurveParams^.TowerParam2.Gamma2.toHexString
else Result:=CurveParams^.TowerParam2.Gamma2.Inverse.toHexString
end;

{*******************************************************************************}
function TKSS36Curve.getH: string;
begin
Result:=CurveParams^.H.toHexString;
end;

{*******************************************************************************}
function TKSS36Curve.getHtw: string;
begin
Result:=CurveParams^.Htw.toHexString;
end;

{*******************************************************************************}
function TKSS36Curve.getLp: string;
begin
Result:=CurveParams^.Lp.toHexString;
end;

{*******************************************************************************}
function TKSS36Curve.getN: string;
begin
Result:=CurveParams^.N.toHexString;
end;

{*******************************************************************************}
function TKSS36Curve.getP: string;
begin
Result:=CurveParams^.P.toHexString;
end;

{*******************************************************************************}
function TKSS36Curve.getR: string;
begin
Result:=CurveParams^.R.toHexString;
end;

function TKSS36Curve.GetRandomG1Point: FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsRandomTorsionPoint;
end;

function TKSS36Curve.GetRandomG2Point: Fp6Point;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsRandomTorsionPoint;
end;

{*******************************************************************************}
function TKSS36Curve.getRtw: String;
begin
Result:=CurveParams^.Rtw.toHexString;
end;

{*******************************************************************************}
function TKSS36Curve.getSigma: String;
begin
Result:=CurveParams^.TowerParam.Sigma.toHexString
end;

{*******************************************************************************}
function TKSS36Curve.getTr: string;
begin
Result:=CurveParams^.Tr.toHexString;
end;

function TKSS36Curve.GetTwm: TTwistModel;
begin
Result:=CurveParams.TwistMode;
end;

{*******************************************************************************}
function TKSS36Curve.getu: String;
begin
Result:=_u.ToHexString;
end;

function TKSS36Curve.HashToG1Point(id: string): FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsTorsionFromString(id);
end;

function TKSS36Curve.HashToG2Point(id: String): Fp6Point;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsTorsionFromString(id);
end;

{*******************************************************************************}
procedure TKSS36Curve.RecomputeParametres;
var tmp:KSS36CurvesParamsDefinition;
    tmph,q,r,power:HLInt;
begin
tmp.u:=_u.toHexString;
tmp.Beta:=beta;
tmp.Sigma:=Sigma;
tmp.Gamma:=Gamma;
tmp.A:=_A.toHexString;
tmp.B:=_B.toHexString;
tmp.SecurityLevel:=SecurityLevel;
ComputeKSS36Parametres(tmp,CurveParams);
end;

{*******************************************************************************}
procedure TKSS36Curve.SetStandardCurveParametres(inParametres:StandardKSS36Curves);
var tmph,q,r,power:HLInt;
begin
case inParametres of
scKSS36at256_1:begin
               ComputeKSS36Parametres(KSS36_256_1_Params,CurveParams);
               _A:=CurveParams.A;
               _B:=_Str_To_LInt(KSS36_256_1_Params.B);
               _u:=_Hex_To_LInt(KSS36_256_1_Params.u);
                end;
scKSS36at256_2:begin
               ComputeKSS36Parametres(KSS36_256_2_Params,CurveParams);
               _A:=CurveParams.A;
               _B:=_Str_To_LInt(KSS36_256_2_Params.B);
               _u:=_Hex_To_LInt(KSS36_256_2_Params.u);
              end;
scKSS36at192_1:begin
               ComputeKSS36Parametres(KSS36_192_1_Params,CurveParams);
               _A:=CurveParams.A;
               _B:=_Str_To_LInt(KSS36_192_1_Params.B);
               _u:=_Hex_To_LInt(KSS36_192_1_Params.u);
              end;
scKSS36at192_2:begin
               ComputeKSS36Parametres(KSS36_192_2_Params,CurveParams);
               _A:=CurveParams.A;
               _B:=_Str_To_LInt(KSS36_192_2_Params.B);
               _u:=_Hex_To_LInt(KSS36_192_2_Params.u);
               end;
end;
VLargeIntegers._Pow_LInt(CurveParams.P,power,12);
VLargeIntegers._Pow_LInt(CurveParams.P,tmph,6);
VLargeIntegers._Sub_LInt(power,tmph,power);
VLargeIntegers._Inc_LInt(power,1);
VLargeIntegers._Div_Mod_LInt(power,CurveParams.R,power,tmph);
q:=power;
Setlength(Decomp,0);
while not VLargeIntegers._IsNull(q) do begin
                        r:=Q mod CurveParams.p;
                        q:=q /CurveParams.p;
                        Setlength(Decomp,length(Decomp)+1);
                        _VHCopy_LInt(r, Decomp[length(Decomp)-1]);
                        end;
if HammingWeight(Lp,true)<HammingWeight(lp,false) then PerferedPoweringMode:=lpmNaf
else PerferedPoweringMode:=lpmBinary;
end;

procedure TKSS36Curve.SetPArams(const Value: StandardKSS36Curves);
var pr:TComponent;
    i:integer;
begin
SetStandardCurveParametres(Value);
pr:=Self.Owner;
for i:=0 to pr.ComponentCount-1 do begin
                                   if (pr.Components[i] is TTreeView1)and (TTreeView1(pr.Components[i]).Curve = self)
                                   then GenerateParamsTreeView(TTreeView(pr.Components[i]));
                                   end;
params:=value;
end;

{*******************************************************************************}
procedure TKSS36Curve.SetA(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TKSS36Curve.SetB(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TKSS36Curve.SetBeta(Value: Integer);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TKSS36Curve.setCoord(value: CoordinatesSystem);
begin
CurveParams^.CoordSys:=value;
end;

{*******************************************************************************}
procedure TKSS36Curve.SetCustomCurveParametres(inParametres: KSS36CurvesParamsDefinition);
begin
ComputeKSS36Parametres(inParametres,CurveParams);
end;

procedure TKSS36Curve.SetGamma(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TKSS36Curve.SetSigma(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TKSS36Curve.Setu(Value: String);
begin
RecomputeParametres;
end;

{********* Pairing function :Compute the pairing according to the specified algorithm***********}
function TKSS36Curve.Paire(P: FpPoint; Q: Fp6Point): Fp36Int;
begin
if (P.Infinity)or(Q.Infinity) then raise Exception.Create('Can''t compute pairing for infinit points ....');
Result:=OptAtePairing(P,Q);
end;

{******************** Define the curve used for the pairing ***********************}
procedure TKSS36Curve.SetCurve(Value: TKSS36Curve);
var tmp,q,r,power:HLInt;
begin
_LoopMode:= PerferedPoweringMode;
SetLoopMode(_LoopMode);
VLargeIntegers._Pow_LInt(CurveParams.P,power,12);
VLargeIntegers._Pow_LInt(CurveParams.P,tmp,6);
VLargeIntegers._Sub_LInt(power,tmp,power);
VLargeIntegers._Inc_LInt(power,1);
VLargeIntegers._Div_Mod_LInt(power,CurveParams.R,power,tmp);
q:=power;
while not VLargeIntegers._IsNull(q) do begin
                        r:=Q mod CurveParams.p;
                        q:=q /CurveParams.p;
                        Setlength(Decomp,length(Decomp)+1);
                        _VHCopy_LInt(r, Decomp[length(Decomp)-1]);
                        end;
end;


{******************** Set the loop mode of the pairing :Binary / Negatve Representation*************}
procedure TKSS36Curve.SetLoopMode(Value: LoopPoweringMode);
begin
_LoopMode:=Value;
if _LoopMode=lpmAuto then _LoopMode:=PerferedPoweringMode;

case _LoopMode of
lpmBinary:begin
          case PairingAlgorithm of
          KSS36OptAte:Loop:=CurveParams.LoopBin;
          end;
          end;
lpmNaf:begin
       case PairingAlgorithm of
       KSS36OptAte:Loop:=CurveParams.LoopNaf;
       end;
       end;
end;
end;

procedure _Pow_FP36_list(const Value: FP36Int; List: LIntArray; var Result: FP36Int);
var  i,j: word;
     tmp,f2,zi: FP36Int;
     nafexpo,t1,t2:LIntArrayForm;
     SizeNaf,Sizebin,max:Integer;
     res:array of FP36Int;
     dec:boolean;
begin
Result.SetTowerParams(Value.Tower);
tmp := Value;
max:=0;
for i:=0 to length(List)-1 do if max<List[i].BitLength then max:=List[i].BitLength;
zi.SetTowerParams(Value.Tower);
Setlength(Res,length(List));
for i:=0 to length(Res)-1 do begin
                             Res[i].SetTowerParams(Value.Tower);
                             if List[i].IsOdd then Res[i]:=value
                             else Res[i].SetToOne;
                             end;
for i :=1 to max-1 do begin
                     //  _Sqr_FP36(tmp, tmp);
                       _Compressed_Sqr_FP36(tmp,tmp);
                       dec:=false;
                       for j:=0 to length(List)-1 do begin
                                                    if _Is_BitSet_At(List[j], i) then begin
                                                                                     if not dec then begin
                                                                                                     _DeCompressed_FP36(tmp,zi);
                                                                                                     dec:=true;
                                                                                                     end;
                                                                                      _Mul_FP36(zi,Res[j], Res[j]);
                                                                                      end;
                                                    end;
                       end;
Result:=Res[0];
for i:=1 to length(res)-1 do begin
                             _Pow_FP36_P_i(Res[i],i,f2);
                             _Mul_FP36(f2 ,result,result);
                             end;
end;
{******************** Final Exponentiation Step :Compute f^((p^36-1)/r) ***********************}
Procedure TKSS36Curve.FinalPowerOptimalAteKSS1(f:Fp36Int;u:LInt;PoweringMode:FpPoweringMode;var Result:Fp36Int);
var f1,f2,f3,exp:Fp36Int;
    i:integer;

begin
for i:=0 to 30 do t[i].SetTowerParams(f.Tower);

///***** Soft Exponentiation ***///
_Inv_FP36(f,f1);
_Conjugate_FP36(f,f2);
_Mul_FP36(f2,f1,f2);
_Pow_FP36_P_i(f2,6,f1);
_Mul_FP36(f2,f1,f1);
///***** Hard Exponentiation ***///
_Pow_FP36_list(f1,Decomp,Result);
end;


{******************** Compute the Optimal Ate Pairing ***********************}
Function TKSS36Curve.OptAtePairing(Pt:FpPoint;Qt:Fp6Point):Fp36Int;
var T,mQt,Q1,Q2:Fp6Point;
    i:integer;
    f,f1,f2:Fp36Int;
    passed:boolean;
begin
f.SetTowerParams(CurveParams.TowerParam2);
f1.SetTowerParams(CurveParams.TowerParam2);
f.SetToOne;
T:=Qt;
Q1:=Qt;
_Neg_Fp6_Point(Qt,mQt);
T.SetPairingPointCoordinates(Pt.X,Pt.Y);
Q1.SetPairingPointCoordinates(Pt.X,Pt.Y);
T.ComputeLigneValue:=true;
Q1.ComputeLigneValue:=true;
passed:=false;
for i:=Length(Loop^)-2 Downto 0 do begin
                                   case CoordinatesSystem of
                                      csAffine:_Double_Affine_Fp6_Point(T,T,true);
                                      csProjective:_Double_Projective_Fp6_Point(T,T,True);
                                   end;
                                   _Sqr_Fp36(f,f);
                                   if not(passed) then begin
                                                       Q1:=T;
                                                       passed:=true;
                                                       end;
                                   //_Sparse_Mul_Fp36(f,T.LineAtP36,f,Curve.TwistMode);
                                   _Mul_Fp36(f,T.LineAtP36,f);
                                   if  Loop^[i]=1 then begin
                                                        case CoordinatesSystem of
                                                          csAffine:_Add_Affine_Fp6_Point(Qt,T,T,True);
                                                          csProjective:_Add_Projective_Fp6_Point(Qt,T,T,True);
                                                          end;
                                                       //_Sparse_Mul_Fp36(f,T.LineAtP36,f,Curve.TwistMode);
                                                       _Mul_Fp36(f,T.LineAtP36,f);
                                                       end
                                   else if (Loop^[i]=-1) then begin
                                                                case CoordinatesSystem of
                                                                  csAffine:_Add_Affine_Fp6_Point(mQt,T,T,True);
                                                                  csProjective:_Add_Projective_Fp6_Point(mQt,T,T,True);
                                                                  end;
                                                              //_Sparse_Mul_Fp36(f,T.LineAtP36,f,Curve.TwistMode);
                                                              _Mul_Fp36(f,T.LineAtP36,f);
                                                              end;
                                   end;
if _IsNeg(CurveParams.u) then begin
                                    _Neg_Fp6_Point(T,T);
                                    _Conjugate_FP36(f,f);
                                    end;
_Pow_FP36_P_i(Q1.LineAtP36,7,f1);
f2:=Q1.LineAtP36;
case CoordinatesSystem of
     csAffine:_Add_Affine_Fp6_Point(Qt,Q1,Q1,True);
     csProjective:_Add_Projective_Fp6_Point(Qt,Q1,Q1,True);
end;
f2:=f2*Q1.LineAtP36;
_Conjugate_FP36(f2,f2);
_Pow_FP36_P_i(f2,1,f2);
_Mul_Fp6_FpPoint(-3*CurveParams.FieldParam.p,Qt,Q1);
case CoordinatesSystem of
      csAffine:_Add_Affine_Fp6_Point(Q1,T,T,True);
      csProjective:_Add_Projective_Fp6_Point(Q1,T,T,True);
    end;
f:=f*f1*f2*T.LineAtP36;
FinalPowerOptimalAteKSS1(f,Qt.CurveParams.u,Fp36PoweringMode,Result);
end;

procedure Register;
begin
  RegisterComponents('Pairings', [TKSS36Curve]);
end;

end.

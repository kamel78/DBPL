unit BLS27Curves;

interface

uses   Vcl.dialogs,Vcl.ComCtrls, System.SysUtils,VLargeIntegers,LargeIntegers,Fp9Arithmetic,Fp3Arithmetic,GeneralTypes,
      System.classes,ECCFp9,ECCFp,Fp27Arithmetic ;


Type
   G2BLS27=Fp9Point;
   GTBLS27=Fp27Int;
   G1BLS27=Fppoint;
   LargeInt=Lint;
   BLS27CurvesPairingAlgos=(BLS27pOptAte);

   BLS27CurvesParamsDefinition=record
                            SecurityLevel:String;
                            u:String;  // the paramater of generation for the BLS27 curve
                            Beta:Integer; // non-square elements of the irredictible polynomial on FP2
                            Sigma:string; // non-square non-cube elements of the irredictible polynomial on FP4
                            Gamma:string; // non-square non-cube elements of the irredictible polynomial on FP8
                            A,B:String; // parametres of the curve   y^2=x^3+Ax+B
                            TwistMode:TTwistModel;
                            GenratorX:String;
                            TwistGeneratorSeed:Word;
                            end;
   ListOfBLS27Params=array of BLS27CurvesParamsDefinition;

 TBLS27Curve=class (TCurve)
                      {Definition of a Curve with respect the the Weistrass model   Y^2=X^3+A*X+B (mod P)}
              private
              _A,_B:Lint;
              _u:Lint;
              _LoopMode:LoopPoweringMode;                     //  Loopmode :Negative/Binary Representation
              _Fp27PoweringMode:FpPoweringMode;             // Final Exponentiation Powering Mode : Beuchat Approach, Karabina Approach, Classical Square and Multiply
              _PairingAlgorithm:BLS27CurvesPairingAlgos;         //  Pairing Algorithme OptAte,R-Ate.....
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
              Procedure FinalPowerOptimalAteBLS1(f:Fp27Int;u:LInt;PoweringMode:FpPoweringMode;var Result:Fp27Int); // Compute f^((p^24-1)/r)
              Function OptAtePairing(Pt:FpPoint;Qt:Fp9Point):Fp27Int;  // Compute Optimal Ate Pairings
              procedure SetLoopMode(Value:LoopPoweringMode);
            public
              SecurityLevel:String;    /// A symbolic identifier of the curve
              PerferedPoweringMode:LoopPoweringMode;
              Loop:PLIntArrayForm;
              property Beta:integer read getBeta write SetBeta;
              property Sigma:string read Getsigma write SetSigma;
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
              procedure SetStandardCurveParametres(inParametres:BLS27CurvesParamsDefinition);
              procedure SetCustomCurveParametres(inParametres:BLS27CurvesParamsDefinition);
              procedure GenerateParamsTreeView(Tree:TTreeView);override;
              function Paire(P:FpPoint;Q:Fp9Point):Fp27Int;
              function GetRandomG1Point:FpPoint;
              function GetRandomG2Point:Fp9Point;
              function GetDefautG1Generator:FpPoint;
              function GetDefautG2Generator:Fp9Point;
              function HashToG1Point(id:string):FpPoint;
              function HashToG2Point(id:String):Fp9Point;
              published
              //property Fp27PoweringMode:FpPoweringMode read _Fp27PoweringMode write _Fp27PoweringMode;
              property PairingAlgorithm:BLS27CurvesPairingAlgos read _PairingAlgorithm write _PairingAlgorithm ;
              property LoopMode:LoopPoweringMode read _LoopMode write SetLoopMode;
              property CoordinatesSystem:CoordinatesSystem read getCoord write SetCoord default csAffine;
            end;

const

  BLS27_256_1_Params:BLS27CurvesParamsDefinition=(SecurityLevel:'192bit';u:'0x1A0000F8';Beta:2;sigma:'0+u*1+u^2*0';A:'0';B:'-2';TwistMode:twMType;GenratorX:'-1';TwistGeneratorSeed:55);

  BLS27ParamsList:array[0..0] of string=('BLS27_256_1_Params');
  BLS27ImplementedPairingAlgos:array[0..0] of string=('Optimal Ate Pairing');

  procedure ComputeBLS27Parametres(Params: BLS27CurvesParamsDefinition; var Result:PtrCurveParams);
  procedure Register;

implementation



procedure ComputeBLS27Parametres(Params: BLS27CurvesParamsDefinition; var Result:PtrCurveParams);
var u9,u6,u3,tmp,tmp1:LInt;
    i:integer;
    m2,_P:HLInt;
    e:Lint;
    gamma,twtmp,twtmp1:Fp9Int;
    tau:array[0..9] of HLint;
       s:string;
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
                Family:=cfBLS27;
                if (Pos('0x',Params.u)<>0)or(Pos('$',Params.u)<>0) then u:=_hex_To_LInt(Params.u)
                else u:=_str_To_LInt(Params.u);
                u9:=u.Sqr.Sqr.Sqr*u;
                r:=(u9.Sqr+u9+1)/3;
                p:=((u-1).Sqr*r)+u;
                Tr:=u+1;
                n:=p+1-tr;
                if (not IsLIntPrime(p)){or(not IsLIntPrime(r)) }then
                    raise Exception.Create('Paramétre t invalide pour la construction de la courbe BLS27...')
                else begin
                     p.Limit:=2*p.Data.i16[-2];
                     new(p.InverseFordivision);
                     p.InverseFordivision^:=(Lint(1) shl (p.Limit*32))/p+1;
                     p.id:='Y';
                     Lp:=(u).Absolute;
                     end;
                A:=_Str_To_LInt(Params.A);
                B:=_Str_To_LInt(Params.B);
                FieldParam:=InitFieldParamsFp3(Params.Beta,p);
                //computing number of points on E(FP9)with cofactor
                Tau[0]:=2;
                Tau[1]:=Tr;
                for i:=1 to 8 do Tau[i+1]:=tr*Tau[i]-p*Tau[i-1];
                //s:=Tau[9].ToDecimalString;
                  _P:=p.Sqr.Sqr.Sqr*p;
                _Sqrt_HLint((4*_P-Tau[9].Sqr)/3,tau[0]);

                m2:=_P+1-((-3*tau[0]-tau[9])/2);
                tau[0]:=m2/r;
                _VHCopy_LInt(tau[0],Htw); // Twist Cofactor
                H:=n/r;


                Rtw:=r;         ///  The true order is m2=r*htw , so we can choose a subgourp of order R (the cofactor is htw) (https://eprint.iacr.org/2005/133.pdf)

                New(TowerParam3);
                TowerParam3^.FieldParam:=FieldParam;
                TowerParam3^.Sigma.SetFieldParams(TowerParam3^.FieldParam);
                TowerParam3.Sigma.SetFormString(Params.Sigma);
                LoopBin^:=Lp.ToIntArray;
                LoopNaf^:=Lp.ToNafArray;
                LoopTateBin^:=n.ToIntArray;
                LoopTateNaf^:=n.ToNafArray;

                TowerParam3.pmod8:=p mod 8;
                BtwFp9.SetTowerParams(TowerParam3);
                AtwFp9.SetTowerParams(TowerParam3);
                gamma.SetTowerParams(TowerParam3);
                gamma.a.SetToZero;
                gamma.b.SetToOne;
                AtwFp9.SetToZero;
                gamma.c.SetToZero;
                TwistMode:=Params.TwistMode;
                if TwistMode=twMType then _Mul_FP_FP9(B,Gamma.Sqr,BtwFp9) else
                _Mul_FP_FP9(B,gamma.Sqr.Inverse,BtwFp9);

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

                // constants for powering Fp9 elements to P
                TowerParam3^.FrobeniusP_Const[0].SetFieldParams(TowerParam3^.FieldParam);
                TowerParam3^.FrobeniusP_Const[0]:=TowerParam3^.Sigma.Pow((P-1)/3);
                TowerParam3^.FrobeniusP_Const[1].SetFieldParams(TowerParam3^.FieldParam);
                TowerParam3^.FrobeniusP_Const[1]:=TowerParam3^.FrobeniusP_Const[0]*TowerParam3^.FrobeniusP_Const[0];
                TowerParam3^.FrobeniusP_Const[2]:=TowerParam3^.Sigma.Pow((P-1)/9);
                TowerParam3^.FrobeniusP_Const[3]:=TowerParam3^.Sigma.Pow((P-1)/9).Sqr;
                FrobeniusMapConstX_Fp9.SetTowerParams(TowerParam3);
                case TwistMode of
                    twDType: begin
                             FrobeniusMapConstX_Fp9:=gamma.Pow((2*p-2)/3);
                             FrobeniusMapConstY_Fp9:=gamma.Pow(p-1);
                             end;
                twMType:begin
                        FrobeniusMapConstX_Fp9:=(gamma.Pow((2*p-2)/3)).Inverse;///
                        FrobeniusMapConstY_Fp9:=(gamma.Pow(p-1)).Inverse;   //////
                        end;
                end;
                                ///   Constructing the base point for G2 Generator
                BLS27TwistGeneratorX.SetTowerParams(TowerParam3);
                BLS27TwistGeneratorY.SetTowerParams(TowerParam3);

                BLS27TwistGeneratorX.a.SetToRandom;
                BLS27TwistGeneratorX.b.SetToRandom;
                BLS27TwistGeneratorX.c.SetToRandom;

                repeat
                _Inc_LInt(BLS27TwistGeneratorX.a.a,1);
                _Sqr_FP9(BLS27TwistGeneratorX,twtmp);
                _Add_FP9(twtmp,AtwFp9,twtmp);
                _Mul_FP9(twtmp,BLS27TwistGeneratorX,twtmp1);
                _Add_FP9(twtmp1,BtwFp9,twtmp);
                until twtmp.IsASquare(BLS27TwistGeneratorY);

                end;
end;

{ TBLS27Curve }

procedure TBLS27Curve.GenerateParamsTreeView(Tree:TTreeView);
var s:string;
    r:real;
    G:FpPoint;Gtw:Fp9Point;
    tmp:TTreeNode;
begin
Tree.Items.Clear;
tmp:=Tree.Items.Add(nil,'Curve');
Tree.Items.AddChild(tmp,'Family : BLS27');
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
Tree.Items.AddChild(tmp,'Equation of the Twist : '+'y^2=x^3+'+Curveparams.BtwFp9.toHexString);
Tree.Items.AddChild(tmp,'Degree of the Twist : Cubic');
if Curveparams.TwistMode=twMType then S:='M-Type' else S:='D-Type';
Tree.Items.AddChild(tmp,'Type of the Twist : '+S);
tmp:=Tree.Items.Add(nil,'Field G2 (Twist) :E''(Fp9)[R]');
Tree.Items.AddChild(tmp,'Order of G2 R = '+CurveParams.Rtw.ToHexString);
Tree.Items.AddChild(tmp,'Size of G2 : '+inttostr(CurveParams.Rtw.BitLength)+' bit');
Tree.Items.AddChild(tmp,'Cofactor of G2 Htw = '+CurveParams.Htw.ToHexString);
tmp:=Tree.Items.Add(nil,'Targted Field (GT):');
Tree.Items.AddChild(tmp,'Extension Field : Fp27');
Tree.Items.AddChild(tmp,'Size of GT : '+inttostr(27*CurveParams.P.BitLength)+' bit');
tmp:=Tree.Items.Add(nil,'Security Level');
Tree.Items.AddChild(tmp,'Security in G1 : '+inttostr(CurveParams.R.BitLength div 2)+' bit').ImageIndex:=2;
Tree.Items.AddChild(tmp,'Security in G2 : '+inttostr(CurveParams.R.BitLength div 2)+' bit').ImageIndex:=2;
r:=1.92* exp((1/3)*ln(27*CurveParams.P.BitLength))*sqr(exp((1/3)*(ln(ln(27*CurveParams.P.BitLength)/ln(2)))));
Tree.Items.AddChild(tmp,'Security in GT (GNFS) : '+inttostr(Round(r))+' bit').ImageIndex:=2;
tmp:=Tree.Items.Add(nil,'Implemented Pairings');
Tree.Items.AddChild(tmp,'Otp-Ate');
tmp:=Tree.Items.Add(nil,'Tower Construction');
if Curveparams.FieldParam.Beta<0 then Tree.Items.AddChild(tmp,'Fp3<u>=ExstensionField<u,|u^3+'+inttostr(Curveparams.FieldParam.Beta)+'>')
else Tree.Items.AddChild(tmp,'Fp3<u>=ExstensionField<u,|u^3-'+inttostr(Curveparams.FieldParam.Beta)+'>');
Tree.Items.AddChild(tmp,'Fp9<v>=ExstensionField<v,|v^3-('+Curveparams.TowerParam3.Sigma.toHexString+')>');
Tree.Items.AddChild(tmp,'Fp27<w>=ExstensionField<t,|w^3-v>');
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


constructor TBLS27Curve.Create;
begin
inherited Create(AOwner);
CurveParams:=nil;
SetStandardCurveParametres(BLS27_256_1_Params);
CoordinatesSystem:=csAffine;
_LoopMode:= PerferedPoweringMode;
SetLoopMode(_LoopMode);
end;

{*******************************************************************************}
destructor TBLS27Curve.Destroy;
begin
Dispose(CurveParams);
inherited;
end;

{*******************************************************************************}
function TBLS27Curve.getA: String;
begin
Result:=_A.ToHexString;
end;

{*******************************************************************************}
function TBLS27Curve.getAtw: String;
begin
Result:=CurveParams^.AtwFp9.toHexString;
end;

{*******************************************************************************}
function TBLS27Curve.getB: String;
begin
Result:=_B.ToHexString;
end;

{*******************************************************************************}
function TBLS27Curve.getBeta:integer;
begin
Result:=CurveParams^.FieldParam.Beta;
end;

{*******************************************************************************}
function TBLS27Curve.getBtw: String;
begin
Result:=CurveParams^.BtwFp9.toHexString;
end;

{*******************************************************************************}
function TBLS27Curve.getCoord: CoordinatesSystem;
begin
Result:=CurveParams^.CoordSys;
end;

function TBLS27Curve.GetDefautG1Generator: FpPoint;
begin
Result.SetToDefaultGenerator;
end;

function TBLS27Curve.GetDefautG2Generator: Fp9Point;
begin
Result.SetToDefaultGenerator;
end;

{*******************************************************************************}
function TBLS27Curve.getH: string;
begin
Result:=CurveParams^.H.toHexString;
end;

{*******************************************************************************}
function TBLS27Curve.getHtw: string;
begin
Result:=CurveParams^.Htw.ToDecimalString;
end;

{*******************************************************************************}
function TBLS27Curve.getLp: string;
begin
Result:=CurveParams^.Lp.toHexString;
end;

{*******************************************************************************}
function TBLS27Curve.getN: string;
begin
Result:=CurveParams^.N.toHexString;
end;

{*******************************************************************************}
function TBLS27Curve.getP: string;
begin
Result:=CurveParams^.P.toHexString;
end;

{*******************************************************************************}
function TBLS27Curve.getR: string;
begin
Result:=CurveParams^.R.toHexString;
end;

{*******************************************************************************}
function TBLS27Curve.GetRandomG1Point: FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsRandomTorsionPoint;
end;

{*******************************************************************************}
function TBLS27Curve.GetRandomG2Point: Fp9Point;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsRandomTorsionPoint;
end;

{*******************************************************************************}
function TBLS27Curve.getRtw: String;
begin
Result:=CurveParams^.Rtw.toHexString;
end;

{*******************************************************************************}
function TBLS27Curve.getSigma: String;
begin
Result:=CurveParams^.TowerParam3.Sigma.toHexString
end;

{*******************************************************************************}
function TBLS27Curve.getTr: string;
begin
Result:=CurveParams^.Tr.toHexString;
end;

{*******************************************************************************}
function TBLS27Curve.GetTwm: TTwistModel;
begin
Result:=CurveParams.TwistMode;
end;

{*******************************************************************************}
function TBLS27Curve.getu: String;
begin
Result:=_u.ToHexString;
end;

{*******************************************************************************}
function TBLS27Curve.HashToG1Point(id: string): FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsTorsionFromString(id);
end;

{*******************************************************************************}
function TBLS27Curve.HashToG2Point(id: String): Fp9Point;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsTorsionFromString(id);
end;

{*******************************************************************************}
procedure TBLS27Curve.RecomputeParametres;
var tmp:BLS27CurvesParamsDefinition;
begin
tmp.u:=_u.toHexString;
tmp.Beta:=beta;
tmp.Sigma:=Sigma;
tmp.A:=_A.toHexString;
tmp.B:=_B.toHexString;
tmp.SecurityLevel:=SecurityLevel;
ComputeBLS27Parametres(tmp,CurveParams);
end;

{*******************************************************************************}
procedure TBLS27Curve.SetStandardCurveParametres(inParametres:BLS27CurvesParamsDefinition);
begin
begin
ComputeBLS27Parametres(inParametres,CurveParams);
     _A:=CurveParams.A;
     _B:=_Str_To_LInt(inParametres.B);
     _u:=_Hex_To_LInt(inParametres.u);
     end;
if HammingWeight(Lp,true)<HammingWeight(lp,false) then PerferedPoweringMode:=lpmNaf
else PerferedPoweringMode:=lpmBinary;
end;

{*******************************************************************************}
procedure TBLS27Curve.SetA(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBLS27Curve.SetB(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBLS27Curve.SetBeta(Value: Integer);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBLS27Curve.setCoord(value: CoordinatesSystem);
begin
if value<>csAffine then begin
                        Messagedlg('Only Affine coordinaes system is implemented for the BLS27 Family',Mtwarning,[mbok],0);
                        CurveParams^.CoordSys:=csAffine;
                        end
else CurveParams^.CoordSys:=value;
end;

{*******************************************************************************}
procedure TBLS27Curve.SetCustomCurveParametres(inParametres: BLS27CurvesParamsDefinition);
begin
ComputeBLS27Parametres(inParametres,CurveParams);
end;

{*******************************************************************************}
procedure TBLS27Curve.SetSigma(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TBLS27Curve.Setu(Value: String);
begin
RecomputeParametres;
end;

{********* Pairing function :Compute the pairing according to the specified algorithm***********}
function TBLS27Curve.Paire(P: FpPoint; Q: Fp9Point): Fp27Int;
begin
if (P.Infinity)or(Q.Infinity) then raise Exception.Create('Can''t compute pairing for infinit points ....');
Result:=OptAtePairing(P,Q);
end;

{******************** Set the loop mode of the pairing :Binary / Negatve Representation*************}
procedure TBLS27Curve.SetLoopMode(Value: LoopPoweringMode);
begin
_LoopMode:=Value;
if _LoopMode=lpmAuto then _LoopMode:=PerferedPoweringMode;

case _LoopMode of
lpmBinary:begin
          case PairingAlgorithm of
          BLS27pOptAte:Loop:=CurveParams.LoopBin;
          end;
          end;
lpmNaf:begin
       case PairingAlgorithm of
       BLS27pOptAte:Loop:=CurveParams.LoopNaf;
       end;
       end;
end;
end;

{******************** Final Exponentiation Step :Compute f^((p^27-1)/r) ***********************}
Procedure TBLS27Curve.FinalPowerOptimalAteBLS1(f:Fp27Int;u:LInt;PoweringMode:FpPoweringMode;var Result:Fp27Int);
var m,m1,m2,m3,tmp0,tmp1:Fp27Int;
    z_1,z3,tmp:Lint;
    tt:TTime;
    h,mm,s,ms:word;
begin
tt:=now;
///***** Soft Exponentiation ***///
_Pow_Fp27_P9(f,m);
_Inv_Fp27(f,m1);
_Mul_Fp27(m1,m,M);
///***** Hard Exponentiation ***///

_Pow_Fp27_P_i(m,8,m1);


_Pow_Fp27(m,u,tmp1);

_Pow_Fp27_P_i(tmp1,7,tmp0);
_Mul_Fp27(tmp0,m1,m1);

_Pow_Fp27(tmp1,u,tmp1);
_Pow_Fp27_P_i(tmp1,6,tmp0);
_Mul_Fp27(tmp0,m1,m1);

_Pow_Fp27(tmp1,u,tmp1);
_Pow_Fp27_P_i(tmp1,5,tmp0);
_Mul_Fp27(tmp0,m1,m1);

_Pow_Fp27(tmp1,u,tmp1);
_Pow_Fp27_P_i(tmp1,4,tmp0);
_Mul_Fp27(tmp0,m1,m1);

_Pow_Fp27(tmp1,u,tmp1);
_Pow_Fp27_P_i(tmp1,3,tmp0);
_Mul_Fp27(tmp0,m1,m1);

_Pow_Fp27(tmp1,u,tmp1);
_Pow_Fp27_P_i(tmp1,2,tmp0);
_Mul_Fp27(tmp0,m1,m1);

_Pow_Fp27(tmp1,u,tmp1);
_Pow_Fp27_P(tmp1,tmp0);
_Mul_Fp27(tmp0,m1,m1);

_Pow_Fp27(tmp1,u,tmp1);
_Mul_Fp27(tmp1,m1,m1);

_HCopy_LInt(u,z_1);
_Dec_LInt(z_1,1);



_Pow_Fp27(m1,z_1,m2);
_Pow_Fp27(m2,z_1,m2);

_Sqr_LInt(u,tmp);
_Mul_LInt(tmp,u,z3);

_Pow_Fp27(m2,z3,m3);
_Pow_Fp27(m3,z3,m3);
_Pow_Fp27(m3,z3,m3);

_Mul_Fp27(m3,m2,m3);
_Pow_Fp27_P9(m2,m2);
_Mul_Fp27(m3,m2,m3);

_Sqr_Fp27(m,Result);
_Mul_Fp27(Result,m,Result);
_Mul_Fp27(Result,m3,Result);
  //tt:=now-tt;
//decodetime(tt,h,mm,s,ms);showmessage(inttostr(h)+inttostr(mm)+inttostr(s)+inttostr(ms));
end;


{******************** Compute the Optimal Ate Pairing ***********************}
Function TBLS27Curve.OptAtePairing(Pt:FpPoint;Qt:Fp9Point):Fp27Int;
var T,mQt:Fp9Point;
    i:integer;
    f,f1,f2:Fp27Int;
    tt:TTime;
    h,m,s,ms:word;
begin
tt:=now;
f.SetTowerParams(CurveParams.TowerParam3);
f.SetToOne;
T:=Qt;
_Neg_Fp9_Point(Qt,mQt);
T.SetPairingPointCoordinates(Pt.X,Pt.Y);
T.ComputeLigneValue:=true;
for i:=Length(Loop^)-2 Downto 0 do begin
                                   case CoordinatesSystem of
                                      csAffine:_Double_Affine_Fp9_Point(T,T,true);
                                      csProjective:;/// Implementing Projective coordinates on BLS27 is not necessary , the gain is is negligible
                                   end;
                                   _Sqr_Fp27(f,f);
                                   _Mul_Fp27(f,T.LineAtP27,f);
                                   if  Loop^[i]=1 then begin
                                                        case CoordinatesSystem of
                                                          csAffine:_Add_Affine_Fp9_Point(Qt,T,T,True);
                                                          csProjective:;/// Implementing Projective coordinates on BLS27 is not necessary , the gain is is negligible
                                                          end;
                                                       _Mul_Fp27(f,T.LineAtP27,f);
                                                       end
                                   else if (Loop^[i]=-1) then begin
                                                                case CoordinatesSystem of
                                                                  csAffine:_Add_Affine_Fp9_Point(T,mqt,T,True);
                                                                  csProjective:;/// Implementing Projective coordinates on BLS27 is not necessary , the gain is is negligible
                                                                  end;
                                                              _Mul_Fp27(f,T.LineAtP27,f);
                                                              end;
                                   end;
tt:=now-tt;
//decodetime(tt,h,m,s,ms);showmessage(inttostr(h)+inttostr(m)+inttostr(s)+inttostr(ms));
FinalPowerOptimalAteBLS1(f,CurveParams.u,pmNormal,Result);
end;

procedure Register;
begin
  RegisterComponents('Pairings', [TBLS27Curve]);
end;


end.

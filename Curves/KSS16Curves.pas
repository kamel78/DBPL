unit KSS16Curves;

interface

uses  Vcl.ComCtrls, System.SysUtils,VLargeIntegers,LargeIntegers,Fp16Arithmetic,Fp8Arithmetic,Fp4Arithmetic
      ,Fp2Arithmetic,GeneralTypes,System.classes,VCL.dialogs,ECCFp4,ECCFp, TreeView1;


Type

   GTKSS16=Fp16Int;
   G1KSS16=FpPoint;
   G2KSS16=Fp4Point;
   LargeInt=Lint;
   KSS16CurvesPairingAlgos=(KSS16pOptAte);
   StandardKSS16Curves=(sc128Kss16_1,sc128Kss16_2_raz,sc128Kss16_3,sc192Kss16_1,sc192Kss16_2);
   KSS16CurvesParamsDefinition=record
                            SecurityLevel:String;
                            u:String;  // the paramater of generation for the KSS16 curve
                            Beta:Integer; // non-square elements of the irredictible polynomial on FP2
                            Sigma:string; // non-square non-cube elements of the irredictible polynomial on FP4
                            Gamma:string; // non-square non-cube elements of the irredictible polynomial on FP8
                            A,B:String; // parametres of the curve   y^2=x^3+Ax+B
                            TwistMode:TTwistModel;
                            GenratorX:String;
                            TwistGeneratorSeed:Word;
                            end;
   ListOfKSS16Params=array of KSS16CurvesParamsDefinition;

 TKSS16CurvePairing=class (TCurve)
                      {Definition of a Curve with respect the the Weistrass model   Y^2=X^3+A*X+B (mod P)}
              private
              _A,_B:Lint;
              _u:Lint;
              _LoopMode:LoopPoweringMode;                     //  Loopmode :Negative/Binary Representation
              _FP16PoweringMode:FpKSS16PoweringMode;             // Final Exponentiation Powering Mode : Beuchat Approach, Karabina Approach, Classical Square and Multiply
              _PairingAlgorithm:KSS16CurvesPairingAlgos;         //  Pairing Algorithme OptAte,R-Ate.....
              Loop:PLIntArrayForm;
              params:StandardKSS16Curves;
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
              Procedure FinalPowerOptimalAteKSS16(ff:Fp16Int;u:LInt;PoweringMode:FpKSS16PoweringMode;var Result:Fp16Int); // Compute f^((p^16-1)/r)
              Function OptAtePairing(Pt:FpPoint;Qt:Fp4Point):Fp16Int;  // Compute Optimal Ate Pairings
              procedure SetLoopMode(Value:LoopPoweringMode);
              procedure SetPArams(const Value: StandardKSS16Curves);
            public
              SecurityLevel:String;    /// A symbolic identifier of the curve
              PerferedPoweringMode:LoopPoweringMode;
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
              procedure SetStandardCurveParametres(inParametres:StandardKSS16Curves);
              procedure SetCustomCurveParametres(inParametres:KSS16CurvesParamsDefinition);
              procedure GenerateParamsTreeView(Tree:TTreeView);override;
              function Paire(P:FpPoint;Q:Fp4Point):Fp16Int;
              function GetRandomG1Point:FpPoint;
              function GetRandomG2Point:Fp4Point;
              function GetDefautG1Generator:FpPoint;
              function GetDefautG2Generator:Fp4Point;
              function HashToG1Point(id:string):FpPoint;
              function HashToG2Point(id:String):Fp4Point;
              published
              property Parametres:StandardKSS16Curves read params write Setparams;
              property FP16PoweringMode:FpKSS16PoweringMode read _FP16PoweringMode write _FP16PoweringMode;
              property PairingAlgorithm:KSS16CurvesPairingAlgos read _PairingAlgorithm write _PairingAlgorithm ;
              property LoopMode:LoopPoweringMode read _LoopMode write SetLoopMode;

            end;

const

 KSS16_128_1_Params:KSS16CurvesParamsDefinition=(SecurityLevel:'128bit';u:'67130955';Beta:2;sigma:'0+u*1';Gamma:'(0+u*0,1+u*0)';A:'2';B:'0';TwistMode:twMType;GenratorX:'-1';TwistGeneratorSeed:55);
 KSS16_128_2_Params:KSS16CurvesParamsDefinition=(SecurityLevel:'128bit';u:'-0x3F87007FF';Beta:-2;sigma:'0+u*1';Gamma:'(0+u*0,1+u*0)';A:'1';B:'0';TwistMode:twMType;GenratorX:'3';TwistGeneratorSeed:55);
 KSS16_Razvan_128_Params:KSS16CurvesParamsDefinition=(SecurityLevel:'128bit';u:'0x6fffc0101';Beta:2;sigma:'0+u*1';Gamma:'(0+u*0,1+u*0)';A:'1';B:'0';TwistMode:twDType;GenratorX:'3';TwistGeneratorSeed:55);
 KSS16_192_1_Params:KSS16CurvesParamsDefinition=(SecurityLevel:'192bit';u:'0x2000004007F7F';Beta:-3;sigma:'0+u*1';Gamma:'(0+u*0,1+u*0)';A:'3';B:'0';TwistMode:twMType;GenratorX:'-1';TwistGeneratorSeed:55);
 KSS16_192_2_Params:KSS16CurvesParamsDefinition=(SecurityLevel:'192bit';u:'-0xDFFFFDFFEFFF';Beta:3;sigma:'0+u*1';Gamma:'(0+u*0,1+u*0)';A:'3';B:'0';TwistMode:twMType;GenratorX:'-1';TwistGeneratorSeed:55);

 KSS16ParamsList:array[0..4] of string=('KSS16_128_1_Params','KSS16_Razvan_128_Params','KSS16_128_2_Params','KSS16_192_1_Params','KSS16_192_2_Params');
 KSS16ImplementedPairingAlgos:array[0..0] of string=('Optimal Ate Pairing');

 procedure ComputeKSS16Parametres(Params: KSS16CurvesParamsDefinition; var Result:PtrCurveParams);
 procedure Register;

implementation



procedure ComputeKSS16Parametres(Params: KSS16CurvesParamsDefinition; var Result:PtrCurveParams);
var u4,u8,u3,tmp,tmp1:LInt;
    i:integer;
    m2:HLInt;
    e:Lint;
    _P,_F,_tau:HLInt;
    twtmp:Fp4Int;
    twtmp1:Fp4Int;
    TwoinFp2:Fp2int;
    tau:array[0..4] of HLint;
    Genrator,Gtmp:Fp4Point;
    BLSTwistBasePointX,BLSTwistBasePointY:Fp4Int;
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
                Family:=cfKSS16;
                if (Pos('0x',Params.u)<>0)or(Pos('$',Params.u)<>0) then u:=_hex_To_LInt(Params.u)
                else u:=_str_To_LInt(Params.u);
                u4:=(u.Sqr).Sqr;
                u8:=u4.Sqr;
                Tr:=(2*u4*u + 41*u + 35)/35;
                r:=(u8 + 48*u4 + 625)/61250;
                p:=(u.Pow(10)+2*u.Pow(9)+5*u.Pow(8)+48*u.Pow(6)+152*u.Pow(5)+240*u.Pow(4)+625*u.Sqr+2398*u+3125)/980;
                n:=p+1-tr;
                h:=n/r;
                if ((not IsLIntPrime(p))or(not IsLIntPrime(r)))or(n mod r <>0) then
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
                //computing number of points on E(FP4)with cofactor
                e:=params.beta;
                if e.IsASqrMod(P,FieldParam.MontgomeryData) then raise Exception.Create(inttostr(Params.beta)+' is a quadratic Residue modulo P......');
                Tau[0]:=2;
                Tau[1]:=tr;
                for i := 1 to 3 do Tau[i+1]:=Tr*Tau[i]-p*Tau[i-1];
                _P:=p.Sqr*p.Sqr;
                _Tau:=Tau[4];
                _F:=(4*_P-_Tau.Sqr);
                _Sqrt_HLint(_F,_F);
                m2:=_P+1+_F;
                _VHCopy_LInt(m2,tmp1);

                Htw:=tmp1/r;    // Twist Cofactor
                Rtw:=r;         ///  The true order is m2=r*htw , so we can choose a subgourp of order R (the cofactor is htw) (https://eprint.iacr.org/2005/133.pdf)

                New(TowerParam);
                TowerParam^.FieldParam:=FieldParam;
                TowerParam^.Sigma.SetFieldParams(TowerParam^.FieldParam);
                TowerParam.Sigma.SetFormString(Params.Sigma);
                     ///   should test if beta and sigma are valide parametres
                TwoinFp2.SetFieldParams(FieldParam);
                _Pow_FP2(TowerParam.Sigma,(Result.P.Sqr-1)/2,TwoinFp2);
                if (TwoinFp2.a=1)and(TwoinFp2.b=0) then raise Exception.Create('Paramétre t invalide pour la construction de la courbe KSS...Sigma is a square');
                {_Pow_FP2(Result.TowerParam.Sigma,(Result.P.Sqr-1)/3,TwoinFp2);
                if (TwoinFp2.a=1)and(TwoinFp2.b=0) then raise Exception.Create('Paramétre t invalide pour la construction de la courbe KSS...Sigma is a Cube');}
                    ///
                TowerParam.pmod8:=p mod 8;
                TwoinFp2.SetFieldParams(TowerParam^.FieldParam);
                TwoinFp2.SetFormString('2+u*0');       // for squaring over FP4
                TowerParam^.inv_2.SetFieldParams(TowerParam^.FieldParam);
                TowerParam^.inv_2:=TwoinFP2.Inverse;

                New(TowerParam2);
                TowerParam2.Tower8paprms:=TowerParam;
                TowerParam2.pmod8:=TowerParam.pmod8;
                TowerParam2.Sigma:=Towerparam.Sigma;
                TowerParam2.FieldParam:=FieldParam;
                TowerParam2.Gamma.SetTowerParams(TowerParam);
                TowerParam2.Gamma.SetFromStrings(Copy(Params.Gamma,2,pos(',',Params.Gamma)-2),Copy(Params.Gamma,pos(',',Params.Gamma)+1,pos(')',Params.Gamma)-pos(',',Params.Gamma)-1));
                BtwFp4.SetTowerParams (TowerParam);
                TwistMode:=Params.TwistMode;

                BtwFp4.SetTowerParams(TowerParam);
                BtwFp4.SetFromStrings('0+u*0','0+u*0');
                AtwFp4.SetTowerParams(TowerParam);
                TwistMode:=Params.TwistMode;
                if TwistMode=twMType then _Mul_FP_FP4(A,TowerParam2.Gamma,AtwFp4) else
                _Mul_FP_FP4(A,TowerParam2.Gamma.Inverse,AtwFp4);

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
                BLSTwistBasePointX.SetTowerParams(TowerParam);
                BLSTwistBasePointY.SetTowerParams(TowerParam);
                BLSTwistBasePointX.a.a.SetToRandom(p,params.TwistGeneratorSeed);
                BLSTwistBasePointX.a.b.SetToRandom(p,params.TwistGeneratorSeed);
                BLSTwistBasePointX.b.a.SetToRandom(p,params.TwistGeneratorSeed);
                BLSTwistBasePointX.b.b.SetToRandom(p,params.TwistGeneratorSeed);

                repeat
                _Inc_LInt(BLSTwistBasePointX.a.a,1);
                _Sqr_FP4(BLSTwistBasePointX,twtmp);
                _Add_FP4(twtmp,AtwFp4,twtmp);
                _Mul_FP4(twtmp,BLSTwistBasePointX,twtmp1);
                _Add_FP4(twtmp1,BtwFp4,twtmp);
                until twtmp.IsASquare;
                _Sqrt_FP4(twtmp,BLSTwistBasePointY);

                TowerParam2^.FrobeniusP_Const_Fp4[0].SetTowerParams(TowerParam);
                TowerParam2^.FrobeniusP_Const_Fp4[0]:=TowerParam2^.Gamma.Pow((P-1)/4);
                for i:=1 to 2 do begin
                                 TowerParam2^.FrobeniusP_Const_Fp4[i].SetTowerParams(TowerParam);
                                 TowerParam2^.FrobeniusP_Const_Fp4[i]:=TowerParam2^.FrobeniusP_Const_Fp4[i-1]*TowerParam2^.FrobeniusP_Const_Fp4[0];
                                 end;
                FrobeniusMapConstX_Fp4.SetTowerParams(TowerParam);
                case TwistMode of
                    twDType: begin
                             FrobeniusMapConstX_Fp4:=TowerParam2^.Gamma.Pow((p-1)/2);
                             FrobeniusMapConstY_Fp4:=TowerParam2^.Gamma.Pow(3*(p-1)/4);
                             end;
                twMType:begin
                        FrobeniusMapConstX_Fp4:=(TowerParam2^.Gamma.Pow((p-1)/2)).Inverse;///
                        FrobeniusMapConstY_Fp4:=(TowerParam2^.Gamma.Pow(3*(p-1)/4)).Inverse;   //////
                        end;
                end;
                TowerParam.FrobeniusPi3xSigmaSqr:=(TowerParam2^.Gamma.Pow((p-3)).a);
                TowerParam.FrobeniusPi3xSigmaSqr:=TowerParam.FrobeniusPi3xSigmaSqr*TowerParam2.Sigma;
                TowerParam2.FrobeniusPi3xSigmaSqr:=TowerParam.FrobeniusPi3xSigmaSqr;
                Genrator.SetCurveParams(Result);
                Gtmp.SetCurveParams(Result);
                Gtmp.X:=BLSTwistBasePointX;
                Gtmp.Y:=BLSTwistBasePointY;
                Gtmp.Z.a.SetFormString('1');
                Gtmp.Z.b.SetFormString('0');
                Gtmp.Infinity:=false;
                _Mul_Fp_Fp4Point(htw,Gtmp,Genrator);
                BLS24TwistGeneratorX:=Genrator.X;
                BLS24TwistGeneratorY:=Genrator.Y;
                end;
end;


{ TKSS16Curve }

procedure TKSS16CurvePairing.GenerateParamsTreeView(Tree:TTreeView);
var s:string;
    r:real;
    G:FpPoint;Gtw:Fp4Point;
    tmp:TTreeNode;
begin
Tree.Items.Clear;
tmp:=Tree.Items.Add(nil,'Curve');
Tree.Items.AddChild(tmp,'Family : KSS16');
if CurveParams.A<0 then  s:='y^2=x^3-'+CurveParams.A.Absolute.ToDecimalString+'*x'
else s:='y^2=x^3+'+CurveParams.A.ToDecimalString+'*x';
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
Tree.Items.AddChild(tmp,'Equation of the Twist : '+'y^2=x^3+'+Curveparams.BtwFp4.toHexString);
Tree.Items.AddChild(tmp,'Degree of the Twist : Quartic');
if Curveparams.TwistMode=twMType then S:='M-Type' else S:='D-Type';
Tree.Items.AddChild(tmp,'Type of the Twist : '+S);
tmp:=Tree.Items.Add(nil,'Field G2 (Twist) :E''(Fp4)[R]');
Tree.Items.AddChild(tmp,'Order of G2 R = '+CurveParams.Rtw.ToHexString);
Tree.Items.AddChild(tmp,'Size of G2 : '+inttostr(CurveParams.Rtw.BitLength)+' bit');
Tree.Items.AddChild(tmp,'Cofactor of G2 Htw = '+CurveParams.Htw.ToHexString);
tmp:=Tree.Items.Add(nil,'Targted Field (GT):');
Tree.Items.AddChild(tmp,'Extension Field : Fp16');
Tree.Items.AddChild(tmp,'Size of GT : '+inttostr(16*CurveParams.P.BitLength)+' bit');
tmp:=Tree.Items.Add(nil,'Security Level');
Tree.Items.AddChild(tmp,'Security in G1 : '+inttostr(CurveParams.R.BitLength div 2)+' bit').ImageIndex:=2;
Tree.Items.AddChild(tmp,'Security in G2 : '+inttostr(CurveParams.R.BitLength div 2)+' bit').ImageIndex:=2;
r:=1.92* exp((1/3)*ln(16*CurveParams.P.BitLength))*sqr(exp((1/3)*(ln(ln(16*CurveParams.P.BitLength)/ln(2)))));
Tree.Items.AddChild(tmp,'Security in GT (GNFS) : '+inttostr(Round(r))+' bit').ImageIndex:=2;
tmp:=Tree.Items.Add(nil,'Implemented Pairings');
Tree.Items.AddChild(tmp,'Otp-Ate');
tmp:=Tree.Items.Add(nil,'Tower Construction');
if Curveparams.FieldParam.Beta<0 then Tree.Items.AddChild(tmp,'Fp2<u>=ExstensionField<u,|u^2+'+inttostr(Curveparams.FieldParam.Beta)+'>')
else Tree.Items.AddChild(tmp,'Fp2<u>=ExstensionField<u,|u^2-'+inttostr(Curveparams.FieldParam.Beta)+'>');
Tree.Items.AddChild(tmp,'Fp4<v>=ExstensionField<v,|v^2-('+Curveparams.TowerParam.Sigma.toHexString+')>');
Tree.Items.AddChild(tmp,'Fp8<w>=ExstensionField<w,|w^2-v>');
Tree.Items.AddChild(tmp,'Fp16<z>=ExstensionField<z,|z^2-w>');
tmp:=Tree.Items.Add(nil,' = '+floattostrf(CurveParams.P.BitLength/CurveParams.R.BitLength,ffGeneral, 2, 4 ));
tmp:=Tree.Items.Add(nil,'Default G1 Generator');
G.SetCurveParams (CurveParams);
G.SetToDefaultGenerator;
Gtw.SetCurveParams(CurveParams);
Gtw.SetToDefaultGenerator;
Tree.Items.AddChild(tmp,G.toHexString);
tmp:=Tree.Items.Add(nil,'Default G2 Generator');
Tree.Items.AddChild(tmp,Gtw.toHexString);
Tree.FullExpand;
Tree.Selected:=Tree.Items[0];
//Tree.SetFocus;
end;

constructor TKSS16CurvePairing.Create;
begin
inherited create(Aowner);
CurveParams:=nil;
SetStandardCurveParametres(sc128Kss16_1);
CoordinatesSystem:=csProjective;
FP16PoweringMode:=pmkssNormal;
SetLoopMode(lpmAuto);
PairingAlgorithm:=KSS16pOptAte;
end;

{*******************************************************************************}
destructor TKSS16CurvePairing.Destroy;
begin
Dispose(CurveParams);
inherited;
end;

{*******************************************************************************}
function TKSS16CurvePairing.getA: String;
begin
Result:=_A.ToHexString;
end;

{*******************************************************************************}
function TKSS16CurvePairing.getAtw: String;
begin
Result:=CurveParams^.AtwFp4.toHexString;
end;

{*******************************************************************************}
function TKSS16CurvePairing.getB: String;
begin
Result:=_B.ToHexString;
end;

{*******************************************************************************}
function TKSS16CurvePairing.getBeta:integer;
begin
Result:=CurveParams^.FieldParam.Beta;
end;

{*******************************************************************************}
function TKSS16CurvePairing.getBtw: String;
begin
Result:=CurveParams^.BtwFp4.toHexString;
end;

{*******************************************************************************}
function TKSS16CurvePairing.getCoord: CoordinatesSystem;
begin
Result:=CurveParams^.CoordSys;
end;

{*******************************************************************************}
function TKSS16CurvePairing.GetDefautG1Generator: FpPoint;
begin
Result.SetToDefaultGenerator;
end;

{*******************************************************************************}
function TKSS16CurvePairing.GetDefautG2Generator: Fp4Point;
begin
Result.SetToDefaultGenerator;
end;

function TKSS16CurvePairing.getGamma: String;
begin
if  TwistMode=twDType then Result:=CurveParams^.TowerParam2.Gamma.toHexString
else Result:=CurveParams^.TowerParam2.Gamma.Inverse.toHexString
end;

{*******************************************************************************}
function TKSS16CurvePairing.getH: string;
begin
Result:=CurveParams^.H.toHexString;
end;

{*******************************************************************************}
function TKSS16CurvePairing.getHtw: string;
begin
Result:=CurveParams^.Htw.toHexString;
end;

{*******************************************************************************}
function TKSS16CurvePairing.getLp: string;
begin
Result:=CurveParams^.Lp.toHexString;
end;

{*******************************************************************************}
function TKSS16CurvePairing.getN: string;
begin
Result:=CurveParams^.N.toHexString;
end;

{*******************************************************************************}
function TKSS16CurvePairing.getP: string;
begin
Result:=CurveParams^.P.toHexString;
end;

{*******************************************************************************}
function TKSS16CurvePairing.getR: string;
begin
Result:=CurveParams^.R.toHexString;
end;

function TKSS16CurvePairing.GetRandomG1Point: FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsRandomTorsionPoint;
end;

function TKSS16CurvePairing.GetRandomG2Point: Fp4Point;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsRandomTorsionPoint;
end;

{*******************************************************************************}
function TKSS16CurvePairing.getRtw: String;
begin
Result:=CurveParams^.Rtw.toHexString;
end;

{*******************************************************************************}
function TKSS16CurvePairing.getSigma: String;
begin
Result:=CurveParams^.TowerParam2.Sigma.toHexString
end;

{*******************************************************************************}
function TKSS16CurvePairing.getTr: string;
begin
Result:=CurveParams^.Tr.toHexString;
end;

function TKSS16CurvePairing.GetTwm: TTwistModel;
begin
Result:=CurveParams.TwistMode;
end;

{*******************************************************************************}
function TKSS16CurvePairing.getu: String;
begin
Result:=_u.ToHexString;
end;

{*******************************************************************************}
function TKSS16CurvePairing.HashToG1Point(id: string): FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsTorsionFromString(id);
end;

{*******************************************************************************}
function TKSS16CurvePairing.HashToG2Point(id: String): Fp4Point;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsTorsionFromString(id);
end;

{*******************************************************************************}
procedure TKSS16CurvePairing.RecomputeParametres;
var tmp:KSS16CurvesParamsDefinition;
begin
tmp.u:=_u.toHexString;
tmp.Beta:=beta;
tmp.Sigma:=Sigma;
tmp.Gamma:=Gamma;
tmp.A:=_A.toHexString;
tmp.B:=_B.toHexString;
tmp.SecurityLevel:=SecurityLevel;
ComputeKSS16Parametres(tmp,CurveParams);
end;

{*******************************************************************************}
procedure TKSS16CurvePairing.SetStandardCurveParametres(inParametres:StandardKSS16Curves);
begin
case inParametres of
sc128Kss16_1:ComputeKSS16Parametres(KSS16_128_1_Params,CurveParams);
sc128Kss16_2_raz:ComputeKSS16Parametres(KSS16_Razvan_128_Params,CurveParams);
sc128Kss16_3:ComputeKSS16Parametres(KSS16_128_2_Params,CurveParams);
sc192Kss16_1:ComputeKSS16Parametres(KSS16_192_1_Params,CurveParams);
sc192Kss16_2:ComputeKSS16Parametres(KSS16_192_2_Params,CurveParams);
end;
_A:=CurveParams.A;
_B:=CurveParams.B;
_u:=CurveParams.u;
if HammingWeight(Lp,true)<HammingWeight(lp,false) then PerferedPoweringMode:=lpmNaf
else PerferedPoweringMode:=lpmBinary;
end;

{*******************************************************************************}
procedure TKSS16CurvePairing.SetA(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TKSS16CurvePairing.SetB(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TKSS16CurvePairing.SetBeta(Value: Integer);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TKSS16CurvePairing.setCoord(value: CoordinatesSystem);
begin
CurveParams^.CoordSys:=value;
end;

{*******************************************************************************}
procedure TKSS16CurvePairing.SetCustomCurveParametres(inParametres: KSS16CurvesParamsDefinition);
begin
ComputeKSS16Parametres(inParametres,CurveParams);
end;

procedure TKSS16CurvePairing.SetGamma(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TKSS16CurvePairing.SetSigma(Value: String);
begin
RecomputeParametres;
end;

{*******************************************************************************}
procedure TKSS16CurvePairing.Setu(Value: String);
begin
RecomputeParametres;
end;

{********* Pairing function :Compute the pairing according to the specified algorithm***********}
function TKSS16CurvePairing.Paire(P: FpPoint; Q: Fp4Point): Fp16Int;
begin
if (P.Infinity)or(Q.Infinity) then raise Exception.Create('Can''t compute pairing for infinit points ....');
case PairingAlgorithm of
  KSS16pOptAte:Result:=OptAtePairing(P,Q);
end;
end;

{******************** Set the loop mode of the pairing :Binary / Negatve Representation*************}
procedure TKSS16CurvePairing.SetLoopMode(Value: LoopPoweringMode);
begin
_LoopMode:=Value;
if _LoopMode=lpmAuto then _LoopMode:=PerferedPoweringMode;
case _LoopMode of
lpmBinary:begin
          case PairingAlgorithm of
          KSS16pOptAte:Loop:=CurveParams.LoopBin;
          end;

          end;
lpmNaf:begin
       case PairingAlgorithm of
       KSS16pOptAte:Loop:=CurveParams.LoopNaf;
       end;
       end;
end;
end;

procedure TKSS16CurvePairing.SetPArams(const Value: StandardKSS16Curves);
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

{******************** Final Exponentiation Step :Compute f^((p^16-1)/r) ***********************}
Procedure TKSS16CurvePairing.FinalPowerOptimalAteKSS16(ff:Fp16Int;u:LInt;PoweringMode:FpKSS16PoweringMode;var Result:Fp16Int); // Compute f^((p^16-1)/r)
var f:Fp16Int;
    t:array[1..14]of Fp16Int;
    tmp,tmp1,tmp2:Fp16Int;
    uplus1:Lint;
begin
_Add_Lint(u,1,uplus1);
///***** Soft Exponentiation ***///
_Conjugate_FP16(ff,t[1]);
_Inv_FP16(ff,t[2]);
_Mul_FP16(t[1],t[2],f);
///***** Hard Exponentiation ***///
///  https://pdfs.semanticscholar.org/e7fe/f0b1af449defb81957e314f71ba54805a131.pdf       (but adapted to reduce the number of variables)
_Sqr_FP16(f,tmp);
_Sqr_FP16(tmp,t[1]);
_Pow_FP16(f,uplus1,tmp);
_Pow_FP16(tmp,uplus1,tmp1);
_Mul_FP16(tmp1,t[1],tmp);
_Pow_FP16(tmp,u,t[2]);
_Sqr_FP16(tmp,tmp2);
_Sqr_FP16(tmp2,tmp2);
_Mul_FP16(tmp2,tmp,tmp2);
_Sqr_FP16(t[1],tmp);
_Sqr_FP16(tmp,tmp);
_Sqr_FP16(tmp,tmp);
_Sqr_FP16(tmp,t[3]);
_Conjugate_FP16(t[1],tmp1);
_Mul_FP16(tmp,tmp1,tmp);
_Sqr_FP16(tmp,tmp);
_Pow_FP16(t[2],u,t[4]);
_Pow_FP16(t[4],u,t[5]);
_Mul_FP16(t[5],tmp,t[6]);
_Pow_FP16(t[6],u,tmp1);
_Sqr_FP16(tmp1,tmp);
_Conjugate_FP16(tmp,tmp);
_Pow_FP16(tmp2,25,t[7]);
_Conjugate_FP16(t[7],t[7]);
_Mul_FP16(tmp,t[7],t[7]);
_Sqr_FP16(t[7],t[8]);
_Sqr_FP16(tmp,tmp);
_Sqr_FP16(tmp,tmp);
_Mul_FP16(tmp,tmp1,tmp);
_Mul_FP16(tmp,t[8],t[8]);
_Pow_FP16(tmp1,u,t[9]);
_Pow_FP16(t[9],u,t[10]);
_Pow_FP16(t[10],u,t[11]);
_Sqr_FP16(t[9],tmp);
_Sqr_FP16(tmp,tmp);
_Pow_FP16(t[2],25,tmp1);
_Sqr_FP16(tmp1,tmp2);
_Mul_FP16(tmp2,tmp1,tmp2);
_Mul_FP16(tmp,tmp2,t[12]);
_Conjugate_FP16(t[12],t[12]);
_Mul_FP16(tmp2,tmp1,tmp2);
_Conjugate_FP16(t[9],t[9]);
_Mul_FP16(tmp,t[9],t[9]);
_Mul_FP16(tmp2,t[9],t[9]);
_Sqr_FP16(t[4],tmp);
_Sqr_FP16(tmp,tmp);
_Mul_FP16(t[4],tmp,tmp);
_Sqr_FP16(tmp,tmp1);
_Mul_FP16(tmp1,t[10],t[13]);
_Sqr_FP16(t[10],tmp2);
_Sqr_FP16(tmp1,t[14]);
_Sqr_FP16(t[14],t[14]);
_Mul_FP16(tmp1,t[14],t[14]);
_Mul_FP16(tmp,t[14],t[14]);
_Mul_FP16(tmp2,t[14],t[14]);
_Pow_FP16(t[5],24,tmp);
_Conjugate_FP16(tmp,tmp);
_Conjugate_FP16(t[11],t[11]);
_Mul_FP16(t[11],tmp,tmp);
_Sqr_FP16(t[3],tmp1);
_Mul_FP16(tmp1,t[3],tmp1);
_Mul_FP16(tmp1,t[1],tmp1);
_Mul_FP16(tmp1,tmp,tmp1);
_Sqr_FP16(t[6],tmp);
_Sqr_FP16(tmp,tmp);
_Mul_FP16(tmp,t[6],tmp);
_Mul_FP16(tmp,t[6],tmp);
_Mul_FP16(tmp,t[6],result);
_Pow_FP16_P(t[12],tmp2);
_Pow_FP16_P_i(tmp1,3,tmp);
_Mul_FP16(tmp,tmp2,tmp2);
_Pow_FP16_P_i(t[9],5,tmp);
_Mul_FP16(tmp,tmp2,tmp2);
_Pow_FP16_P_i(Result,7,tmp);
_Mul_FP16(tmp,tmp2,Result);
_Pow_FP16_P_i(t[7],2,tmp);
_Pow_FP16_P_i(t[8],6,tmp2);
_Mul_FP16(tmp,tmp2,tmp);
_Mul_FP16(tmp,Result,Result);
_Pow_FP16_P_i(t[13],4,tmp);
_Mul_FP16(tmp,Result,Result);
_Mul_FP16(t[14],Result,Result);
end;

{******************** Compute the Optimal Ate Pairing ***********************}
Function TKSS16CurvePairing.OptAtePairing(Pt:FpPoint;Qt:Fp4Point):Fp16Int;
var T,Q1,Q2,mQt:Fp4Point;
    i:integer;
    f,lQQ,f2,f3:Fp16Int;
    Passed:boolean;
    aa:lint;
begin
f.SetTowerParams(CurveParams.TowerParam2);
f.SetToOne;
T:=Qt;
_Neg_Fp4_Point(Qt,mQt);
T.SetPairingPointCoordinates(Pt.X,Pt.Y);
Q1.SetPairingPointCoordinates(Pt.X,Pt.Y);
Q2.SetPairingPointCoordinates(Pt.X,Pt.Y);
T. ComputeLigneValue:=true;
Passed:=false;
for i:=Length(Loop^)-2 Downto 0 do begin
                        case CoordinatesSystem of
                          csAffine:_Double_Affine_Fp4_Point(T,T,true);
                          csProjective:_Double_projective_Fp4_Point(T,T,True);
                        end;
                        if not passed then begin
                                           lQQ:=T.LineAtP16^;
                                           passed:=true;
                                           end;
                        _Sqr_FP16(f,f);
                        _Mul_FP16(f,T.LineAtP16^,f);
                        if  Loop^[i]=1 then begin
                                         case CoordinatesSystem  of
                                         csAffine:_Add_Affine_Fp4_Point(Qt,T,T,True);
                                         csProjective:_Add_Projective_Fp4_Point(Qt,T,T,True);
                                         end;
                                         _Mul_FP16(f,T.LineAtP16^,f);
                                         end
                        else if (Loop^[i]=-1) then begin
                                                case CoordinatesSystem of
                                                csAffine:_Add_Affine_Fp4_Point(mQt,T,T,True);
                                                csProjective:_Add_Projective_Fp4_Point(mQt,T,T,True);
                                                end;
                                                _Mul_FP16(f,T.LineAtP16^,f);
                                                end;
                        end;
_Frobenius_Map_i(Qt,1,Q1);
if _IsNeg(CurveParams.u) then begin
                                    _Neg_Fp4_Point(T,T);
                                    _Conjugate_FP16(f,f);
                                    end;

case CoordinatesSystem of
    csAffine:_Add_Affine_Fp4_Point(Q1,T,T,True);
    csProjective:_Add_Projective_Fp4_Point(Q1,T,T,True);
end;
_Mul_FP16(f,T.LineAtP16^,f);
_Pow_FP16_P_i(f,3,f3);
_Mul_FP16(f3,lQQ,f);
FinalPowerOptimalAteKSS16(f,Qt.CurveParams.u,FP16PoweringMode,Result);
end;

procedure Register;
begin
  RegisterComponents('Pairings', [TKSS16CurvePairing]);
end;

end.

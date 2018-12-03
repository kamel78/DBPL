unit SuperSingularCurves;

interface

uses Vcl.ComCtrls, LargeIntegers,Fp12Arithmetic, Fp2Arithmetic,System.SysUtils, GeneralTypes, ECCFp2,ECCFp,System.Classes, TreeView1;

Type
      StandardSSCurves=(scSS90,scSS66,scSS40);
      SSCurvesPairingAlgos=(ssTate);
      LargeInt=Lint;
      GSS=FpPoint;
      GtSS=Fp2Int;
      SSCurvesParamsDefinition=record
                             Identifier:String;
                             A,B,   // parametres of the curve   y^2=x^3+Ax+B
                             P,    //  the prime Modulo
                             N,     // Total number of the points #E(Fp)
                             R,    // Order of the sub-group on the curve (number of the pointes)
                             H:String;    // Cofactor of the curve H=P div R
                             Beta:Integer;
                             end;

       TSuperSingularCurvePairing=Class (TCurve)
                           private
                           _CoordSystem:CoordinatesSystem;
                           _LoopMode:LoopPoweringMode;                     //  Loopmode :Negative/Binary Representation
                           _PairingAlgorithm:SSCurvesPairingAlgos;         //  Pairing Algorithme Tate,Weil.....
                           Loop:PLIntArrayForm;
                           params:StandardSSCurves;

                           procedure setCoord(value:CoordinatesSystem);
                           procedure SetLoopMode(Value:LoopPoweringMode);
                           procedure FinalPowerTateFp2(x:Fp2Int;h:LInt; var result:Fp2Int);
                           Function  TateParing(Pt,Qt:FpPoint):Fp2Int;

                           function getA:String;
                           function getB:String;
                           function getP:String;
                           function getH:String;
                           function getN:String;
                           function getR:String;
                           procedure SetPArams(const Value: StandardSSCurves);
                           public
                           Identifier:String;
                           property CoordinatesSystem:CoordinatesSystem read _CoordSystem write SetCoord;
                           property A:string read getA;
                           property B:string read getB;
                           property P:string read GetP;
                           property H:string read GetH;
                           property N:string read getN;
                           property R:String read getR;
                           constructor Create(AOwner : TComponent); override;
                           destructor destroy;
                           procedure SetStandardCurveParametres(inParametres:StandardSSCurves);
                           procedure SetCustomCurveParametres(inParametres:SSCurvesParamsDefinition);
                           procedure GenerateRandomCurveParametres(pOrderSize,rOrderSize:integer);
                           procedure GenerateParamsTreeView(Tree:TTreeView);Override;
                           function Paire(P:FpPoint;Q:FpPoint):Fp2Int;
                           function GetRandomGPoint:FpPoint;
                           function GetDefautGGenerator:FpPoint;
                           function HashToG1Point(id:string):FpPoint;
                           published
                           property PairingAlgorithm:SSCurvesPairingAlgos read _PairingAlgorithm write _PairingAlgorithm ;
                           property LoopMode:LoopPoweringMode read _LoopMode write SetLoopMode;
                           property Parametres:StandardSSCurves read params write Setparams;

                           End;
Const
      SUPER_SINGULAR_60bit:SSCurvesParamsDefinition=(Identifier:'SS-512';A:'1';B:'0'; P:'0x400000040000000000000000800000080000000000000000000000000000008080000800000000000000000100000000000000000000000000000000000000FF';
                                                      N:'0x40000004000000000000000080000008000000000000000000000000000000808000080000000000000000010000000000000000000000000000000000000100';
                                                      R:'0x8000000800000000000000000000000000000000000000000000000000000001';
                                                      H:'0x8000000000000000000000010000000000000000000000000000000000000100';
                                                      Beta:-1);

      SUPER_SINGULAR_40bit:SSCurvesParamsDefinition=(Identifier:'SS-168';A:'1';B:'0';P:'$000000920000000000000000000000000000000002480123';
                                                      N:'$000000920000000000000000000000000000000002480124';
                                                      R:'$8000000000000000000000000000000000020001';
                                                      H:'$00000124';
                                                      Beta:-1);
      SUPER_SINGULAR_90bit:SSCurvesParamsDefinition=(Identifier:'SS-1024';A:'1';B:'0'; P:'$40000000004000000000002000000000200000000000000000000000000000008000000000800000000000000000000000000000'+'0000000000000000000000000000000000000000000000000000020000000000000000000001000000000000000000000000000000000000000003FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF';
                                                      N:'$400000000040000000000020000000002000000000000000000000000000000080000000008000000000000000000000000000000000000000000000000000000000000000000'+'0000000000000000200000000000000000000010000000000000000000000000000000000000000040000000000000000000000000000000000';
                                                      R:'$8000000000000000000000400000000000000000000000000000000000000001';
                                                      H:'$800000000080000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000040000000000000000000000000000000000';
                                                      Beta:-1);

   SSParamsList:array[0..2] of string=('SUPER_SINGULAR_40bit','SUPER_SINGULAR_60bit','SUPER_SINGULAR_60bit');
   SSImplementedPairingAlgos:array[0..0] of string=('Tate');

    Function SSParamsToStringForm(p:CurveParams):String;
    procedure Register;

implementation

Function SSParamsToStringForm(p:CurveParams):String;
begin
result:='(A:'''+p.A.ToHexString+''';'+#13+#10+'B:'''+P.B.ToHexString+''';'+#13+#10+
        'P:'''+'$'+p.P.ToHexString+''';'+#13+#10+'N:'''+'$'+p.N.ToHexString+''';'+#13+#10+
        'R:'''+'$'+p.R.ToHexString+''';'+#13+#10+'H:'''+'$'+p.H.ToHexString+''')';
end;

{ TSuperSingularCurve }

constructor TSuperSingularCurvePairing.Create(AOwner : TComponent);
begin
inherited Create(AOwner);
SetStandardCurveParametres(scSS90);
_PairingAlgorithm:=ssTate;
_LoopMode:=lpmBinary;
new(Loop);
end;

destructor TSuperSingularCurvePairing.destroy;
begin
Dispose(CurveParams);
inherited;
end;

procedure TSuperSingularCurvePairing.GenerateRandomCurveParametres(pOrderSize,
  rOrderSize: integer);
var xR,xP,xH: LInt;
    alpha,beta,hOrderSize:integer;
    tmp:SSCurvesParamsDefinition;
begin
hOrderSize:=pOrderSize-rOrderSize;
xR:=1;
xH:=1;
xR:=xR shl (rOrderSize-1)+1;
xH:=xH shl (hOrderSize-1);
Randomize;
alpha:=1;
repeat
  _Set_LInt_BitAt(xR,alpha,false);
  alpha:=random(rOrderSize-1);
  _Set_LInt_BitAt(xR,alpha,true);
  until IsLIntPrime(xR);
alpha:=1;
beta:=1;
repeat
 _Set_LInt_BitAt(xH,alpha,false);
 _Set_LInt_BitAt(xH,beta,false);
 alpha:=random(hOrderSize-2);
 beta:=random(hOrderSize-2);
 _Set_LInt_BitAt(xH,alpha,true);
 _Set_LInt_BitAt(xH,beta,true);
 xP:=(xR*xH)-1;
 until IsLIntPrime(xP) and(xP mod 4=3);
tmp.A:='1';
tmp.B:='0';
tmp.P:=xP.ToHexString;
tmp.R:=xR.ToHexString;
tmp.H:=xH.ToHexString;
tmp.N:=(xP+1).ToHexString;
SetCustomCurveParametres(tmp);
end;

function TSuperSingularCurvePairing.getA: String;
begin
Result:=CurveParams.A.ToHexString
end;

function TSuperSingularCurvePairing.getB: String;
begin
Result:=CurveParams.B.ToHexString
end;

function TSuperSingularCurvePairing.GetDefautGGenerator: FpPoint;
begin
Result:=HashToG1Point('default');
end;

function TSuperSingularCurvePairing.getH: String;
begin
Result:=CurveParams.H.ToHexString
end;

function TSuperSingularCurvePairing.getN: String;
begin
Result:=CurveParams.N.ToHexString
end;

function TSuperSingularCurvePairing.getP: String;
begin
Result:=CurveParams.P.ToHexString
end;


procedure TSuperSingularCurvePairing.GenerateParamsTreeView(Tree:TTreeView);
var s:string;
    r:real;
    tmp:TTreeNode;
begin
Tree.Items.Clear;
tmp:=Tree.Items.Add(nil,'Curve');
Tree.Items.AddChild(tmp,'Family : Super-Singular');
if CurveParams.A=1 then s:='y^2=x^3+x'
else s:='y^2='+CurveParams.A.ToDecimalString+'*x^3+x';
if not(_IsNull(CurveParams.B)) then s:=s+Curveparams.B.ToDecimalString;
Tree.Items.AddChild(tmp,'Equation : '+s);
tmp:=Tree.Items.Add(nil,'Field Fp');
Tree.Items.AddChild(tmp,'Prime P = '+CurveParams.P.ToHexString );
Tree.Items.AddChild(tmp,'Size of Fp : '+inttostr(CurveParams.P.BitLength)+' bit');
tmp:=Tree.Items.Add(nil,'Field E(Fp)');
Tree.Items.AddChild(tmp,'Order (#E(Fp)) N = '+CurveParams.N.ToHexString );
tmp:=Tree.Items.Add(nil,'Base Fields G : E(Fp)[R]');
Tree.Items.AddChild(tmp,'Order of G R = '+CurveParams.R.ToHexString);
Tree.Items.AddChild(tmp,'Size of G : '+inttostr(CurveParams.R.BitLength)+' bit');
Tree.Items.AddChild(tmp,'Cofactor of G H = '+CurveParams.H.ToHexString);
tmp:=Tree.Items.Add(nil,'Targted Field (GT):');
Tree.Items.AddChild(tmp,'Extension Field : Fp2');
Tree.Items.AddChild(tmp,'Size of GT : '+inttostr(2*CurveParams.P.BitLength)+' bit');
tmp:=Tree.Items.Add(nil,'Security Level');
Tree.Items.AddChild(tmp,'Security in G : '+inttostr(CurveParams.R.BitLength div 2)+' bit').ImageIndex:=2;
r:=1.92* exp((1/3)*ln(2*CurveParams.P.BitLength))*sqr(exp((1/3)*(ln(ln(2*CurveParams.P.BitLength)/ln(2)))));
Tree.Items.AddChild(tmp,'Security in GT (GNFS) : '+inttostr(Round(r))+' bit').ImageIndex:=2;
tmp:=Tree.Items.Add(nil,'Implemented Pairings');
Tree.Items.AddChild(tmp,'Tate');
Tree.FullExpand;
Tree.Selected:=Tree.Items[0];
//Tree.SetFocus;
end;

function TSuperSingularCurvePairing.getR: String;
begin
Result:=CurveParams.R.ToHexString
end;

function TSuperSingularCurvePairing.GetRandomGPoint: FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsRandomTorsionPoint;
end;


function TSuperSingularCurvePairing.HashToG1Point(id: string): FpPoint;
begin
Result.SetCurveParams(CurveParams);
Result.SetAsTorsionFromString(id);
end;

procedure TSuperSingularCurvePairing.setCoord(value: CoordinatesSystem);
begin
_CoordSystem:=value;
if CurveParams<>nil then CurveParams.CoordSys:=value;

end;

procedure TSuperSingularCurvePairing.SetCustomCurveParametres(inParametres: SSCurvesParamsDefinition);
begin
if CurveParams=nil then begin
                        new(CurveParams);
                        new(CurveParams.LoopBin);
                        new(CurveParams.LoopNaf);
                        end;
With CurveParams^ do begin
                     Family:=cfSuperSingular;
                     Identifier:=inParametres.Identifier;
                     A:=inParametres.A;
                     B:=inParametres.B;
                     P:=inParametres.P;
                     R:=inParametres.R;
                     H:=inParametres.H;
                     N:=inParametres.N;
                     New(FieldParam);
                     FieldParam^.Beta:=inParametres.Beta;
                     FieldParam^.p:=P;
                     FieldParam^.p_1_div_2:=(p-1) shr 1;
                     FieldParam^.inv_2:=LInt(2).InversModulo(p);
                     New(FieldParam.MontgomeryData);
                     InitMontgomeryStruct(P,FieldParam.MontgomeryData);
                     TowerParam:=nil;
                     CoordSys:=CoordinatesSystem;
                     LoopNaf^:=R.ToNafArray;
                     LoopBin^:=R.ToIntArray;
                    end;
end;

procedure TSuperSingularCurvePairing.SetStandardCurveParametres(inParametres: StandardSSCurves);
begin
case inParametres of
scSS90:SetCustomCurveParametres(SUPER_SINGULAR_90bit);
scSS66:SetCustomCurveParametres(SUPER_SINGULAR_60bit);
scSS40:SetCustomCurveParametres(SUPER_SINGULAR_40bit);
end;
end;

procedure TSuperSingularCurvePairing.FinalPowerTateFp2(x:Fp2Int;h:LInt; var result:Fp2Int);
var t:array[0..3] of LInt;
begin
Result.Field:=x.Field;
_Pow_FP2(x,h,Result);                     //Result:=x.Pow(h);
_Sqr_LInt(Result.a,t[0]);                 //t[0]:=Result.a.Sqr;
_Sqr_LInt(Result.b,t[1]);                 //t[1]:=Result.b.Sqr;
_Add_LInt(t[0],t[1],t[2]);                //t[2]:=((t[1]+t[2]) mod x.Field.p).InversModulo(x.Field.p);
_Mod_LInt(t[2],x.Field.p,t[2]);
_Inv_Mod_LInt(t[2],x.Field.p,t[2]);
_Mul_LInt(Result.a,Result.b,t[3]);        //Result.b:=x.Field.p-((2*Result.a*Result.b*t[3]) mod x.Field.p);(instructions al altered for optimization)
_Mul_LInt(t[3],t[2],Result.b);
_Shl_LInt(Result.b,1);
_Sub_LInt(x.Field.p,Result.b,Result.b);
_Mod_LInt(Result.b,x.Field.p,Result.b);
_Sub_LInt(t[0],t[1],Result.a);            // Result.a:=((t[1]-t[2])*t[3]) mod x.Field.p; (instructions al altered for optimization)
_Mul_LInt(Result.a,t[2],t[0]);
_Mod_LInt(t[0],x.Field.p,Result.a);
end;


Function  TSuperSingularCurvePairing.TateParing(Pt,Qt:FpPoint):Fp2Int;
var Z:FpPoint;
    i:Word;
    Lpq:Fp2Int;
    DistQt:Fp2Point;
    minusPt:FpPoint;
begin
Result.SetFieldParams(CurveParams.FieldParam);
Result.a:=1;
Result.b:=0;
Lpq.SetFieldParams(CurveParams.FieldParam);
DistQt.SetCurveParams(CurveParams);
DistQt.X.a:=Curveparams.P-Qt.X;
DistQt.X.b:=0;
DistQt.Y.a:=0;
DistQt.Y.b:=Qt.Y;
minusPt:=-pt;
pt.Lambda.Data.i32[-1]:=0;
pt.C.Data.i32[-1]:=0;
z.ComputeLineAtQPFp12:=false;
z.ComputeLineAtQPFp6:=false;
pt.ComputeLineAtQPFp12:=false;
pt.ComputeLineAtQPFp6:=false;
if not((Pt.Infinity)or(Qt.Infinity)) then
      begin
      Z:=Pt;
      for i:=length(Loop^)-2 downto 0 do begin
                                        Result:=Result.Sqr;
                                        case CoordinatesSystem of
                                              csAffine:_Double_Affine_Fp_Point(Z,Z,true);
                                              csJacobian:_Double_Jacobian_Fp_Point(Z,Z,true);
                                              csProjective:_Double_Projective_Fp_Point(Z,Z,true);
                                              end;

                                        if not _IsNull(Z.Lambda) then Lpq:=(DistQt.Y-Z.Lambda*DistQT.X)-Z.C
                                        else Lpq:=DistQt.X-Z.C;
                                        Result:=Result*Lpq;
                                        if Loop^[i]<>0 then begin
                                                           if Loop^[i]=1 then begin
                                                                             case CoordinatesSystem of
                                                                             csAffine:_Add_Affine_Fp_Point(pt,Z,Z,true);
                                                                             csJacobian:_Add_Jacobian_Fp_Point(pt,Z,Z,true);
                                                                             csProjective:_Add_Projective_Fp_Point(Z,pt,Z,true);
                                                                             end;

                                                                             end
                                                           else begin
                                                                case CoordinatesSystem of
                                                                csAffine:_Add_Affine_Fp_Point(minusPt,Z,Z,true);
                                                                csJacobian:begin
                                                                           _Add_Jacobian_Fp_Point(minusPt,Z,Z,true);
                                                                           _Jacobian_To_Affine_FpPoint(z,z);
                                                                           end;
                                                                csProjective:_Add_Projective_Fp_Point(Z,minuspt,Z,true);
                                                                end;
                                                                end;
                                                           if not _IsNull(Z.Lambda) then Lpq:=(DistQt.Y-Z.Lambda*DistQt.X)-Z.C
                                                           else Lpq:=DistQt.X-Z.C;
                                                           Result:=Result*Lpq;
                                                           end;
                                        end;
     FinalPowerTateFp2(Result,curveparams.H,Result);
     end;
end;

function TSuperSingularCurvePairing.Paire(P: FpPoint; Q: FpPoint): Fp2Int;
begin
if (P.Infinity)or(Q.Infinity) then raise Exception.Create('Can''t compute pairing for infinit points ....');
case PairingAlgorithm of
ssTate:Result:=TateParing(P,Q);
end;

end;

procedure TSuperSingularCurvePairing.SetLoopMode(Value: LoopPoweringMode);
begin
_LoopMode:=Value;
case value of
lpmBinary:Loop:=CurveParams.LoopBin;
lpmNaf:Loop:=CurveParams.LoopNaf;
end;
end;


procedure TSuperSingularCurvePairing.SetPArams(const Value: StandardSSCurves);
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

procedure Register;
begin
  RegisterComponents('Pairings', [TSuperSingularCurvePairing]);
end;

end.

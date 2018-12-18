unit ECCFp9;

interface

uses Fp9Arithmetic, System.SysUtils, LargeIntegers, Fp27Arithmetic,Fp3Arithmetic, HashFunctions, GeneralTypes, BNCurves,
     vcl.dialogs;
type
{$M 204857600,2048576000}
  Fp9Point = record
    // Definition of a Point on EC Over Fp9
  private
    CurveParams: PtrCurveParams;
  public
    X, Y, Z: Fp9Int;
    ComputeLigneValue: boolean;
    CurrentXp, CurrentYp: LInt;   // Used to store Coordinates of the E(Fp)'s point during Pairings
    LineAtP27: Fp27Int;             // Value of the line between the txo pairing points during Miller loop (On Fp27)
    Infinity: boolean;
    procedure SetCurveParams(Value: PtrCurveParams; ComputLmd: boolean = false);
    procedure SetPointValue(Value: Fp9Int);
    function ToDecimalString: String;
    function ToHexString: String;
    procedure SetAsRandomPoint;
    procedure SetAsRandomTorsionPoint;
    procedure SetToDefaultGenerator;
    procedure SetAsTorsionFromHash(hash: TBytes);
    procedure SetAsTorsionFromString(s:String);
    procedure SetPairingPointCoordinates(PtX, PtY: LInt);
    function CompressToArray: TBytes;
    function IsOnTheCurve: boolean;
    Procedure DeCompressFromArray(a: TBytes);
    function FrobeniusMap(pow:integer): Fp9Point;
    class operator Add(Left, Right: Fp9Point): Fp9Point;
    class operator Subtract(Left, Right: Fp9Point): Fp9Point;
    class operator Multiply(Left: LInt; Right: Fp9Point): Fp9Point;
    class operator Equal(const Left, Right: Fp9Point): boolean;
    class operator NotEqual(const Left, Right: Fp9Point): boolean;
    class operator Negative(const Value: Fp9Point): Fp9Point;
  end;

  procedure _Add_Affine_Fp9_Point(Left, Right: Fp9Point; var Result: Fp9Point;  ComputeLambda: boolean = false);
  procedure _Double_Affine_Fp9_Point(const Value: Fp9Point; var Result: Fp9Point;  ComputeLambda: boolean = false);
  procedure _Double_Fp9_Point(Value: Fp9Point; var Result: Fp9Point;  CoordSys: CoordinatesSystem = csAffine);
  procedure _Add_Fp9_Point(Left, Right: Fp9Point; var Result: Fp9Point;  CoordSys: CoordinatesSystem = csAffine);
  procedure _Sub_Fp9_Point(Left, Right: Fp9Point; var Result: Fp9Point;  CoordSys: CoordinatesSystem = csAffine);
  procedure _Evaluate_Fp9_Point(Value: Fp9Int; var Result: Fp9Point);
  procedure _Neg_Fp9_Point(Value: Fp9Point; var Result: Fp9Point);
  function _Is_OnCurve_Fp9Point(Value: Fp9Point): boolean;
  procedure _Frobenius_Map_i(const Value: Fp9Point; power:integer;var Result: Fp9Point);
  procedure _Mul_Fp9_FpPoint(const Left: LInt; const Right: Fp9Point;  var Result: Fp9Point; CoordSys: CoordinatesSystem = csAffine);
  function _Are_Equals_Fp9Points(Left, Right: Fp9Point): boolean;


implementation

{ ******************************************************************************* }
            /// Procedures for Elliptic Curves Arithmetic over Fp9
{ ********************************************************************-*********** }

    { **********   Add two Fp9 Points  using Affine coordinates ***************** }
procedure _Add_Affine_Fp9_Point(Left, Right: Fp9Point; var Result: Fp9Point;ComputeLambda: boolean = false);
var  tmp,A,B,C,D,E,F,t3:Fp9Int;
begin
if Left.Infinity then Result := Right
else if Right.Infinity then Result := Left
else begin
     if _Equals_FP9(Left.X, Right.X) then begin
                                          if _Equals_Fp9(Left.Y, Right.Y) then _Double_Affine_Fp9_Point(Left, Result)
                                          else begin
                                               Result.Z.a.SetToZero;
                                               Result.Z.b.SetToZero;
                                               Result.Z.c.SetToZero;
                                               Result.SetCurveParams(Right.CurveParams, ComputeLambda);
                                               Result.Infinity := true;
                                               end;
                                          end
     else begin
          Result.SetCurveParams(Left.CurveParams, ComputeLambda);
          _Sub_FP9(Right.X,Left.X,A);
          _Inv_FP9(A,A);
          _Sub_FP9(Right.Y,Left.Y,B);
          _Mul_FP9(B,A,B);
          _Sqr_FP9(B,tmp);
          _Sub_FP9(tmp,Right.X,tmp);
          _Sub_FP9(tmp,Left.X,tmp);
          _Neg_FP9(Right.Y,Result.Y);
          _Sub_FP9(Right.X,tmp,t3);
          _Mul_FP9(t3,B,t3);
          _Sub_FP9(t3,Right.Y,Result.Y);
          Result.X:=tmp;
          Result.Z.a.SetToZero;
          Result.Z.b.SetToZero;
          Result.Z.c.SetToZero;
          if Result.ComputeLigneValue then begin
                                           _Sqr_FP9(Result.X,t3);
                                           _Mul_FP9(B,Result.Y,C);
                                           _Add_FP9(t3,C,D);
                                           _Mul_FP_FP9(Right.CurrentYp,B,E);

                                           if Left.CurveParams.TwistMode=twDType then begin

                                                                                      end
                                            else begin
                                                 _Mul_FP_FP9(Right.CurrentXp,Result.X,Result.LineAtP27.c);
                                                 Result.LineAtP27.b.SetToZero;
                                                 _Sqr_LInt(Right.CurrentXp,Result.LineAtP27.b.b.a);
                                                 _Mod_LInt(Result.LineAtP27.b.b.a,Right.CurveParams.FieldParam.p,Result.LineAtP27.b.b.a);
                                                 _Mul_Fp9_By_V(E,E);
                                                 _Sub_FP9(D,E,Result.LineAtP27.a);
                                                 end;
          end;
          end;
     end;
end;

      { **********   Double an Fp9 Point using Affine coordinates ***************** }
procedure _Double_Affine_Fp9_Point(const Value: Fp9Point; var Result: Fp9Point;ComputeLambda: boolean = false);
var t1,t2,t3,tmp, a, b, c, D, E, F,G,H: Fp9Int;
begin
if (Value.Infinity) or (Value.Y.IsZero) or (Value.Z.IsZero) then Result.Infinity := true
else begin
     Result.Infinity := false;
     Result.SetCurveParams(Value.CurveParams, ComputeLambda);
     _Sqr_Fp9(Value.X, t1);
     _Add_FP9(t1,t1,A);
     _Add_FP9(A,t1,A);
     _Add_FP9(Value.Y,Value.Y,B);
     _Inv_FP9(B,C);
     _Mul_FP9(A,C,D);
     _Sqr_FP9(D,tmp);
     _Sub_FP9(tmp,Value.X,tmp);
     _Sub_FP9(tmp,Value.X,tmp);
     _Neg_FP9(value.y,result.Y);
     _Sub_FP9(value.X,tmp,t1);
     _Mul_FP9(t1,D,t1);
     _Add_FP9(t1,Result.Y,Result.Y);
     result.X:=tmp;
     if Result.ComputeLigneValue then begin
                                      _Sqr_FP9(tmp,t3);
                                      _Mul_FP9(Result.Y,D,E);
                                      _Add_FP9(t3,E,F);
                                      _Mul_FP_FP9(value.CurrentYp, D,G);
                                      _Mul_FP_FP9(value.CurrentXp,Result.X,H);
                                      if Value.CurveParams.TwistMode=twDType then begin

                                      end
                                      else begin
                                           Result.LineAtP27.c:=H;
                                           Result.LineAtP27.b.SetToZero;
                                           _Sqr_LInt(Value.CurrentXp,Result.LineAtP27.b.b.a);
                                           _Mod_LInt(Result.LineAtP27.b.b.a,Value.CurveParams.FieldParam.p,Result.LineAtP27.b.b.a);
                                           _Mul_Fp9_By_V(G,G);
                                           _Sub_FP9(F,G,Result.LineAtP27.a);
                                           end;
                                      end;
  end;
end;



      { **********   Add two Fp9 Points ***************** }
procedure _Add_Fp9_Point(Left, Right: Fp9Point; var Result: Fp9Point;CoordSys: CoordinatesSystem = csAffine);
begin
case CoordSys of
    csAffine:_Add_Affine_Fp9_Point(Left, Right, Result);
    csProjective:;/// Implementing Projective coordinates on BLS27 is not necessary , the gain is is negligible
end;
end;

      { **********   Double an Fp9 Points ***************** }
procedure _Double_Fp9_Point(Value: Fp9Point; var Result: Fp9Point;
  CoordSys: CoordinatesSystem = csAffine);
begin
case CoordSys of
  csAffine:_Double_Affine_Fp9_Point(Value, Result);
  csProjective:;/// Implementing Projective coordinates on BLS27 is not necessary , the gain is is negligible
end;
end;

      { **********   Substract two Fp9 Points ***************** }
procedure _Sub_Fp9_Point(Left, Right: Fp9Point; var Result: Fp9Point;CoordSys: CoordinatesSystem = csAffine);
var  tmp: Fp9Point;
begin
_Neg_Fp9_Point(Right, tmp);
case CoordSys of
    csAffine:_Add_Affine_Fp9_Point(Left, tmp, Result);
    csProjective:;/// Implementing Projective coordinates on BLS27 is not necessary , the gain is is negligible
end;
end;

      { **********   Convert from Projective to Affine for an Fp9 Point ***** }
      // (X/Z,Y/Z)
procedure _Projective_To_Affine_Fp9Point(Value: Fp9Point; var Result: Fp9Point);
var  t1: Fp9Int;
begin
if (Value.Infinity) or (Value.Z.IsZero) then Result.Infinity := true
else begin
     t1 := Value.Z.Inverse;
     Result.X := (Value.X * t1);
     Result.Y := (Value.Y * t1);
     Result.Z.SetToZero;
     Result.Z.a.a:=1;
     end;
end;
      { **********   Convert from Projective to Affine for an Fp9 Point ***** }
      // (X/Z,Y/Z^2)
procedure _Projective2_To_Affine_Fp9Point(Value: Fp9Point; var Result: Fp9Point);
var  t1, t2: Fp9Int;
begin
if (Value.Infinity) or (Value.Z.IsZero) then Result.Infinity := true
else begin
     t1 := Value.Z.Inverse;
     t2 := t1.Sqr;
     Result.X := (Value.X * t1);
     Result.Y := (Value.Y * t2);
     Result.Z.SetToZero;
     Result.Z.a.a:=1;
     end;
end;


    { **********   Negation of an Fp9 Point ************* }
procedure _Neg_Fp9_Point(Value: Fp9Point; var Result: Fp9Point);
begin
Result := Value;
_Neg_Fp9(Result.Y, Result.Y);
end;


    { **********   Multiply a scalar with an Fp9 Point ******* }
procedure _Mul_Fp9_FpPoint(const Left: LInt; const Right: Fp9Point;var Result: Fp9Point; CoordSys: CoordinatesSystem = csAffine);
// Right Should be different than Result
var  i: integer;
begin
Result.CurveParams := Right.CurveParams;
if  Right.CurveParams <>nil then Result.CurveParams:=Right.CurveParams;
if (Right.Infinity) or (_IsNull(Left)) then Result.Infinity := true
  else begin
       Result.Infinity := true;
       if Left = 2 then begin
                        case CoordSys of
                         csAffine:_Double_Affine_Fp9_Point(Right, Result);
                         csProjective:;/// Implementing Projective coordinates on BLS27 is not necessary , the gain is is negligible
                        end;
                        end
       else for i := Left.BitLength - 1 downto 0 do begin
                                                    case CoordSys of
                                                    csAffine:_Double_Affine_Fp9_Point(Result,Result);
                                                    csProjective:;/// Implementing Projective coordinates on BLS27 is not necessary , the gain is is negligible
                                                    end;
                                                    if _Is_BitSet_At(Left,i) then begin
                                                                                  case CoordSys of
                                                                                    csAffine:_Add_Affine_Fp9_Point(Right, Result, Result);
                                                                                    csProjective:;/// Implementing Projective coordinates on BLS27 is not necessary , the gain is is negligible
                                                                                  end;
                                                                                  end;
                                                    end;
       if CoordSys=csProjective then _Projective_To_Affine_Fp9Point(Result, Result);
       if _IsNeg(Left) then _Neg_Fp9_Point(Result,Result);
       end;
end;

      { **********   Multiply a scalar with an Fp9 Point (Negative Representation of the Exponent)******* }
procedure _Mul_NAF_Fp9_FpPoint(const Left: LInt; const Right: Fp9Point;var Result: Fp9Point; CoordSys: CoordinatesSystem = csAffine);
var  i: integer;
     loop: LIntArrayForm;
     NegLeft: Fp9Point;
begin
Result.CurveParams := Right.CurveParams;
_Neg_Fp9_Point(Right, NegLeft);
if (Right.Infinity) or (_IsNull(Left)) then Result.Infinity := true
else  begin
      Result.Infinity := true;
      if Left = 2 then _Double_Fp9_Point(Result, Result, CoordSys)
      else begin
           loop := Left.ToNafArray;
           for i := Length(loop) - 1 downto 0 do begin
                                                 _Double_Fp9_Point(Result, Result, CoordSys);
                                                 if loop[i] = 1 then _Add_Fp9_Point(Right, Result, Result, CoordSys);
                                                 if loop[i] = -1 then _Add_Fp9_Point(NegLeft, Result, Result, CoordSys);
                                                 end;
           _Projective_To_Affine_Fp9Point(Result, Result);
           end;
      end;
end;

      { **********   Compute Frobenius Map for an Fp9 Point ******* }
procedure _Frobenius_Map_i(const Value: Fp9Point; power:integer;var Result: Fp9Point);
var i:integer;
begin
if Value.Infinity then Result.Infinity:=true
else begin
     Result.SetCurveParams(Value.CurveParams);
     Result.X:=Value.X;
     Result.Y:=Value.Y;
     Result.Z:=Value.Z;
     for i:=0 to power-1 do begin
                          _Mul_Fp9(Result.CurveParams.FrobeniusMapConstX_Fp9,Result.X.toPowerP,Result.X);
                          _Mul_Fp9(Result.CurveParams.FrobeniusMapConstY_Fp9,Result.Y.toPowerP,Result.Y);
                          end;
    end;
end;
      { **********   Test if an Fp Point is on the Curve *********** }
function _Is_OnCurve_Fp9Point(Value: Fp9Point): boolean;
var tmp, tmp1: Fp9Int;
begin
if Value.Infinity then Result := true
else begin
     _Sqr_Fp9(Value.X, tmp);
     _Add_Fp9(tmp, Value.CurveParams.AtwFp9, tmp);
     _Mul_Fp9(tmp, Value.X, tmp1);
     _Add_Fp9(tmp1, Value.CurveParams.BtwFp9, tmp);
     _Sqr_Fp9(Value.Y, tmp1);
     Result:=(tmp1=tmp)and (tmp.IsASquare(tmp));
     end;
end;

      { **********   Test if two Fp9 Points are equal ******* }
function _Are_Equals_Fp9Points(Left, Right: Fp9Point): boolean;
begin
Result := (_Equals_Fp9(Left.X, Right.X)) and (_Equals_Fp9(Left.Y, Right.Y));
end;

      { **********   Evaluate an Fp9 Point from X value ****** }
procedure _Evaluate_Fp9_Point(Value: Fp9Int; var Result: Fp9Point);
var  tmp, tmp1: Fp9Int;
begin
Result.X := Value;
Result.Infinity := false;
_Sqr_Fp9(Result.X, tmp);
_Add_Fp9(tmp, Result.CurveParams.AtwFp9, tmp);
_Mul_Fp9(tmp, Result.X, tmp1);
_Add_Fp9(tmp1, Result.CurveParams.BtwFp9, tmp);
if tmp.IsASquare(Result.y) then  begin
                       Result.Z.a.a := 1;
                       Result.Z.a.b := 0;
                       Result.Z.b.a := 0;
                       Result.Z.b.b := 0;
                       Result.Z.c.a := 0;
                       Result.Z.c.b := 0;
                       end
else Result.Infinity := true;
end;

{ ******************************************************************************* }
          // Definitions of an Fp9 Point operators and functions
{ ******************************************************************************* }

{ ******************************************************************************* }
procedure Fp9Point.SetCurveParams(Value: PtrCurveParams;
  ComputLmd: boolean = false);
begin
CurveParams := Value;
if CurveParams.TowerParam3 <> nil then LineAtP27.SetTowerParams(Value.TowerParam3);
X.SetTowerParams(Value.TowerParam3);
Y.SetTowerParams(Value.TowerParam3);
Z.SetTowerParams(Value.TowerParam3);
ComputeLigneValue := ComputLmd;
Infinity := false;
end;

{ ******************************************************************************* }
function Fp9Point.FrobeniusMap(pow:integer): Fp9Point;
begin
_Frobenius_Map_i(Self,pow, Result);
end;

{ ******************************************************************************* }
function Fp9Point.IsOnTheCurve: boolean;
begin
Result := _Is_OnCurve_Fp9Point(Self)
end;

{ ******************************************************************************* }
procedure Fp9Point.SetAsRandomPoint;
begin
Infinity := false;
Randomize;
repeat
    X.a.SetToRandom;
    X.b.SetToRandom;
    X.c.SetToRandom;
    SetPointValue(X);
until not Infinity;
end;

{ ******************************************************************************* }
procedure Fp9Point.SetAsTorsionFromHash(hash: TBytes);
var i, j: integer;
    h1, h2: TBytes;
    tmp: Fp9Point;
begin
j:=0;
CurrentXp:=0;
CurrentYp:=0;
tmp.SetCurveParams(Self.CurveParams);
with tmp do
repeat
  h1:=SHA256BytesHash(hash);
  h1[0]:=j;
  h2:=SHA256BytesHash(h1);
  Infinity:=false;
  X.a.a.Data.i32[length(h1) div 4]:=0;
  for i:=0 to length(h1)-1 do X.a.a.Data.i8[i]:=h1[i];
  X.a.a.Data.i16[-1]:=0;
  X.a.a.Data.i16[-2]:=(length(h1) div 4)+1;
  X.a.a:=X.a.a mod CurveParams.p;
  X.a.b.Data.i32[length(h1) div 4]:=0;
  for i:=0 to length(h1)-1 do X.a.b.Data.i8[i]:=h2[i];
  X.a.b.Data.i16[-1]:=0;
  X.a.b.Data.i16[-2]:=(length(h1) div 4)+1;
  X.a.b:=X.a.b mod CurveParams.p;
  X.a.c:=0;
  X.b.a:=0;
  X.b.b:=0;
  X.b.c:=0;
  X.c.a:=0;
  X.c.b:=0;
  X.c.c:=0;
  SetPointValue(X);
  j:=j+1;
until not Infinity;
_Mul_Fp9_FpPoint(Self.CurveParams.Htw, tmp,Self);
Self.Z.SetToZero;
Self.z.a.a:=1;
end;

{ ******************************************************************************* }
procedure Fp9Point.SetAsTorsionFromString(s:String);
begin
Self.SetAsTorsionFromHash(SHA256StringHash(s));
end;

{ ******************************************************************************* }
procedure Fp9Point.SetAsRandomTorsionPoint;
var tmp: Fp9Point;
begin
tmp.SetCurveParams(CurveParams);
Infinity := false;
Randomize;
{X.a.SetToZero;
X.a.a:=1;
X.b.SetToZero;
X.b.a:=1;}
repeat
//X.b.a:=3;
X.a.SetToRandom;
X.b.SetToRandom;
X.c.SetToRandom;
SetPointValue(X);
if not Infinity then begin
                     //_Fast_Mul_Htw_Fp9Point(Self,tmp);
                     _Mul_Fp9_FpPoint(CurveParams.Htw,Self,tmp);
                     Self := tmp;
                     Self.Z.SetToZero;
                     Self.z.a.a:=1;
                     end;
until not Infinity;
{_Mul_Fp9_FpPoint(CurveParams.Rtw, Self, tmp);   // For Testing Torsion point Validity
if not tmp.Infinity then showmessage('problem! G2');}
end;

{ ******************************************************************************* }
procedure Fp9Point.SetPairingPointCoordinates(PtX, PtY: LInt);
begin
CurrentXp := PtX;
CurrentYp := PtY;
end;

{ ******************************************************************************* }
procedure Fp9Point.SetPointValue(Value: Fp9Int);
begin
_Evaluate_Fp9_Point(Value, Self);
end;

{ ******************************************************************************* }
procedure Fp9Point.SetToDefaultGenerator;
begin
Self.X:=CurveParams.BLS27TwistGeneratorX;
Self.Y:=CurveParams.BLS27TwistGeneratorY;
Self.Z.a.SetFormString('1');
Self.Z.b.SetFormString('0');
Self.Infinity:=false;
end;

{ ******************************************************************************* }
function Fp9Point.ToHexString: String;
begin
if Infinity then Result := 'Infinity'
else Result := '(' + X.ToHexString + ',' + Y.ToHexString + ')';
end;

{ ******************************************************************************* }
function Fp9Point.ToDecimalString: String;
begin
if Infinity then Result := 'Infinity'
else Result := '(' + X.ToDecimalString + ',' + Y.ToDecimalString + ')';
end;

{ ******************************************************************************* }
class operator Fp9Point.Add(Left, Right: Fp9Point): Fp9Point;
begin
_Add_Affine_Fp9_Point(Left, Right, Result);
end;

{ ******************************************************************************* }
class operator Fp9Point.Multiply(Left: LInt; Right { Q } : Fp9Point): Fp9Point;
begin
_Mul_Fp9_FpPoint(Left, Right, Result);
end;

{ ******************************************************************************* }
class operator Fp9Point.Subtract(Left, Right: Fp9Point): Fp9Point;
begin
_Sub_Fp9_Point(Left, Right, Result);
end;

{ ******************************************************************************* }
class operator Fp9Point.Negative(const Value: Fp9Point): Fp9Point;
begin
_Neg_Fp9_Point(Value, Result);
end;

{ ******************************************************************************* }
function Fp9Point.CompressToArray: TBytes;
var L:integer;
begin
if Infinity then begin
                 Setlength(result,0);
                 exit;
                 end;
L:=X.a.Field.p.Data.i32[-1]*4*9;
Setlength(Result,L);
Move(X.a.a.Data.i8[0],Result[0],Length(Result) div 9);
Move(X.a.b.Data.i8[0],Result[L div 9],Length(Result) div 9);
Move(X.a.c.Data.i8[0],Result[L div 9 *2],Length(Result)div 9);
Move(X.b.a.Data.i8[0],Result[L div 9 *3],Length(Result) div 9);
Move(X.b.b.Data.i8[0],Result[L div 9 *4],Length(Result) div 9);
Move(X.b.c.Data.i8[0],Result[L div 9 *5],Length(Result) div 9);
Move(X.c.a.Data.i8[0],Result[L div 9 *6],Length(Result) div 9);
Move(X.c.b.Data.i8[0],Result[L div 9 *7],Length(Result) div 9);
Move(X.c.c.Data.i8[0],Result[L div 9 *8],Length(Result) div 9);
end;

{ ******************************************************************************* }
procedure Fp9Point.DeCompressFromArray(a: TBytes);
var L:integer;
begin
L:=length(a) div 9;
Move(a[0],X.a.a.Data.i8[0],L);
Move(a[L],X.a.b.Data.i8[0],L);
Move(a[L*2],X.a.c.Data.i8[0],L);
Move(a[L*3],X.b.a.Data.i8[0],L);
Move(a[L*4],X.b.b.Data.i8[0],L);
Move(a[L*5],X.b.c.Data.i8[0],L);
Move(a[L*6],X.c.a.Data.i8[0],L);
Move(a[L*7],X.c.b.Data.i8[0],L);
Move(a[L*8],X.c.c.Data.i8[0],L);
X.a.a.Data.i16[-1]:=0;
X.a.b.Data.i16[-1]:=0;
X.a.c.Data.i16[-1]:=0;
X.b.a.Data.i16[-1]:=0;
X.b.b.Data.i16[-1]:=0;
X.b.c.Data.i16[-1]:=0;
X.c.a.Data.i16[-1]:=0;
X.c.b.Data.i16[-1]:=0;
X.c.c.Data.i16[-1]:=0;
X.a.a.Data.i16[-2]:=L div 4;
X.a.b.Data.i16[-2]:=L div 4;
X.a.c.Data.i16[-2]:=L div 4;
X.b.a.Data.i16[-2]:=L div 4;
X.b.b.Data.i16[-2]:=L div 4;
X.b.c.Data.i16[-2]:=L div 4;
X.c.a.Data.i16[-2]:=L div 4;
X.c.b.Data.i16[-2]:=L div 4;
X.c.c.Data.i16[-2]:=L div 4;
SetPointValue(X);
end;

class operator Fp9Point.Equal(const Left, Right: Fp9Point): boolean;
begin
Result := _Are_Equals_Fp9Points(Left, Right);
end;

{ ******************************************************************************* }
class operator Fp9Point.NotEqual(const Left, Right: Fp9Point): boolean;
begin
Result := not _Are_Equals_Fp9Points(Left, Right);
end;

end.


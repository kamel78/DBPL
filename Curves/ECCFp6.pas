unit ECCFp6;

interface

uses Fp6Arithmetic, System.SysUtils, LargeIntegers, Fp36Arithmetic,Fp3Arithmetic, HashFunctions, GeneralTypes, BNCurves,
     vcl.dialogs;
type
{$M 204857600,2048576000}
  Fp6Point = record
    // Definition of a Point on EC Over Fp6
    public
    CurveParams: PtrCurveParams;
    X, Y, Z: Fp6Int;
    ComputeLigneValue: boolean;
    CurrentXp, CurrentYp: LInt;   // Used to store Coordinates of the E(Fp)'s point during Pairings
    LineAtP36: Fp36Int;             // Value of the line between the txo pairing points during Miller loop (On Fp36)
    Infinity: boolean;
    procedure SetCurveParams(Value: PtrCurveParams; ComputLmd: boolean = false);
    procedure SetPointValue(Value: Fp6Int);
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
    function FrobeniusMap(pow:integer): Fp6Point;
    class operator Add(Left, Right: Fp6Point): Fp6Point;
    class operator Subtract(Left, Right: Fp6Point): Fp6Point;
    class operator Multiply(Left: LInt; Right: Fp6Point): Fp6Point;
    class operator Equal(const Left, Right: Fp6Point): boolean;
    class operator NotEqual(const Left, Right: Fp6Point): boolean;
    class operator Negative(const Value: Fp6Point): Fp6Point;
  end;

  procedure _Add_Projective_Fp6_Point(Left { Q } , Right { T } : Fp6Point;  var Result: Fp6Point; ComputeLambda: boolean = false);
  procedure _Double_Projective_Fp6_Point(const Value: Fp6Point;  var Result: Fp6Point; ComputeLambda: boolean = false);
  procedure _Add_Affine_Fp6_Point(Left, Right: Fp6Point; var Result: Fp6Point;  ComputeLambda: boolean = false);
  procedure _Double_Affine_Fp6_Point(const Value: Fp6Point; var Result: Fp6Point;  ComputeLambda: boolean = false);
  procedure _Double_Fp6_Point(Value: Fp6Point; var Result: Fp6Point;  CoordSys: CoordinatesSystem = csAffine);
  procedure _Add_Fp6_Point(Left, Right: Fp6Point; var Result: Fp6Point;  CoordSys: CoordinatesSystem = csAffine);
  procedure _Sub_Fp6_Point(Left, Right: Fp6Point; var Result: Fp6Point;  CoordSys: CoordinatesSystem = csAffine);
  procedure _Evaluate_Fp6_Point(Value: Fp6Int; var Result: Fp6Point);
  procedure _Neg_Fp6_Point(Value: Fp6Point; var Result: Fp6Point);
  function _Is_OnCurve_Fp6Point(Value: Fp6Point): boolean;
  procedure _Frobenius_Map_i(const Value: Fp6Point; power:integer;var Result: Fp6Point);
  procedure _Mul_Fp6_FpPoint(const Left: LInt; const Right: Fp6Point;  var Result: Fp6Point; CoordSys: CoordinatesSystem = csAffine);
  procedure _Projective_To_Affine_Fp6Point(Value: Fp6Point; var Result: Fp6Point);  //(X/Z,Y/Z)
  function _Are_Equals_Fp6Points(Left, Right: Fp6Point): boolean;


implementation

{ ******************************************************************************* }
            /// Procedures for Elliptic Curves Arithmetic over Fp6
{ ********************************************************************-*********** }
var  xA,xB,xC,t0,t1,Qx0_,Qx1,Qx1_,Qx2,Qx2_,Qx3,Qx3_,
     Qx4,Qx4_,Qx5,Qx5_,Qx6,Qx6_,Qx7,Qx7_,Qx8,Qx8_: Fp6Point;
    { **********   Add two Fp6 Points  using Affine coordinates ***************** }
procedure _Add_Affine_Fp6_Point(Left, Right: Fp6Point; var Result: Fp6Point;ComputeLambda: boolean = false);
var  t: array [0 .. 4] of Fp6Int;
     Valx,Valy:Fp6Int;
begin
if Left.Infinity then Result := Right
else if Right.Infinity then Result := Left
else begin
     if _Equals_Fp6(Left.X, Right.X) then begin
                                          if _Equals_Fp6(Left.Y, Right.Y) then _Double_Affine_Fp6_Point(Left, Result)
                                          else begin
                                               Result.Z.a.a := 1;
                                               Result.Z.a.b := 0;
                                               Result.Z.b.a := 0;
                                               Result.Z.b.b := 0;
                                               Result.Z.c.a:=0;
                                               Result.Z.c.b:=0;
                                               Result.SetCurveParams(Right.CurveParams, ComputeLambda);
                                               Result.Infinity := true;
                                               end;
                                          end
     else begin
          Result.SetCurveParams(Left.CurveParams, ComputeLambda);
          ValX:=Right.X;
          ValY:=Right.Y;
          _Sub_Fp6(Right.Y, Left.Y, t[0]);
          _Sub_Fp6(Right.X, Left.X, t[1]);
          _Inv_Fp6(t[1], t[2]);
          _Mul_Fp6(t[0], t[2], t[4]); // Lambda
          _Sqr_Fp6(t[4], t[3]);
          _Sub_Fp6(t[3], Left.X, Result.X);
          _Sub_Fp6(Result.X, Right.X, Result.X);
          _Sub_Fp6(Left.X, Result.X, t[2]);
          _Mul_Fp6(t[4], t[2], Result.Y);
          _Sub_Fp6(Result.Y, Left.Y, Result.Y);
          Result.Z.a.a := 1;
          Result.Z.a.b := 0;
          Result.Z.b.a := 0;
          Result.Z.b.b := 0;
          Result.Z.c.a:=0;
          Result.Z.c.b:=0;
          if Result.ComputeLigneValue then begin
          if Left.CurveParams.TwistMode=twDType then begin
                                                     Result.LineAtP36.a.a.SetToZero;
                                                     _Sub_LInt(Right.CurveParams.FieldParam.p,Right.CurrentYp,Result.LineAtP36.a.a.a.a);
                                                     _Mul_Fp6(t[4],Valx,t[3]);
                                                     _Sub_Fp6(Valy,t[3],Result.LineAtP36.a.b);
                                                     _Mul_FP_Fp6(Right.CurrentXp,t[4],Result.LineAtP36.b.a);
                                                     Result.LineAtP36.b.b.SetToZero;
                                                     Result.LineAtP36.c.a.SetToZero;
                                                     Result.LineAtP36.c.b.SetToZero;
                                                     end
          else begin
               Result.LineAtP36.a.b.SetToZero;
               _Sub_LInt(Right.CurveParams.FieldParam.p,Right.CurrentYp,Result.LineAtP36.a.b.a.a);
               _Mul_Fp6(t[4],Valx,t[3]);
               _Sub_Fp6(Valy,t[3],Result.LineAtP36.a.a);
               _Mul_FP_Fp6(Right.CurrentXp,t[4],Result.LineAtP36.c.a);
               Result.LineAtP36.b.b.SetToZero;
               Result.LineAtP36.b.a.SetToZero;
               Result.LineAtP36.c.b.SetToZero;
               end;
          end;
          end;
     end;
end;

      { **********   Double an Fp6 Point using Affine coordinates ***************** }
procedure _Double_Affine_Fp6_Point(const Value: Fp6Point; var Result: Fp6Point;ComputeLambda: boolean = false);
var t: array [0 .. 3] of Fp6Int;
    ValX, ValY: Fp6Int;
begin
if (Value.Infinity) or (Value.Y.IsZero) then Result.Infinity := true
else  begin
      Result.Infinity := false;
      Result.SetCurveParams(Value.CurveParams, ComputeLambda);
      ValX := Value.X;
      ValY := Value.Y;
      _Sqr_Fp6(Value.X, t[0]);
      _Mul_FP_Fp6(LInt(3), t[0], t[1]);
       _Add_Fp6(t[1],Value.CurveParams.AtwFp6 ,t[1]);   /// Atw is Zero on BLS curves
      _Mul_FP_Fp6(LInt(2), Value.Y, t[3]);
      _Inv_Fp6(t[3], t[3]);
      _Mul_Fp6(t[1], t[3], t[2]); // Lambda
      _Sqr_Fp6(t[2], t[3]);
      _Sub_Fp6(t[3], Value.X, t[3]);
      _Sub_Fp6(t[3], Value.X, Result.X);
      _Sub_Fp6(ValX, Result.X, t[3]);
      _Mul_Fp6(t[2], t[3], Result.Y);
      _Sub_Fp6(Result.Y, ValY, Result.Y);
      if Result.ComputeLigneValue then begin
      if Value.CurveParams.TwistMode=twDType then begin
                                                  Result.LineAtP36.a.a.SetToZero;
                                                  _Sub_LInt(Value.CurveParams.FieldParam.p,Value.CurrentYp,Result.LineAtP36.a.a.a.a);
                                                  _Mul_Fp6(t[2],Valx,t[3]);
                                                  _Sub_Fp6(Valy,t[3],Result.LineAtP36.a.b);
                                                  _Mul_FP_Fp6(Value.CurrentXp,t[2],Result.LineAtP36.b.a);
                                                  Result.LineAtP36.b.b.SetToZero;
                                                  Result.LineAtP36.c.a.SetToZero;
                                                  Result.LineAtP36.c.b.SetToZero;
                                                  end
      else begin
           Result.LineAtP36.a.b.SetToZero;
           _Sub_LInt(Value.CurveParams.FieldParam.p,Value.CurrentYp,Result.LineAtP36.a.b.a.a);
           _Mul_Fp6(t[2],Valx,t[3]);
           _Sub_Fp6(Valy,t[3],Result.LineAtP36.a.a);
           _Mul_FP_Fp6(Value.CurrentXp,t[2],Result.LineAtP36.c.a);
           Result.LineAtP36.b.b.SetToZero;
           Result.LineAtP36.b.a.SetToZero;
           Result.LineAtP36.c.b.SetToZero;
           end;
      end;
      end;
end;

      { **********   Double an Fp6 Point (Projective coordinates (X/Z,Y/Z))***************** }
      // Costello et al's   https://cryptojedi.org/papers/edate-20100614.pdf
procedure _Double_Projective_Fp6_Point(const Value: Fp6Point;var Result: Fp6Point; ComputeLambda: boolean = false);
var tmp,L00,L01,L10, a, b, c, D, E, F, G: Fp6Int;
begin
if (Value.Infinity) or (Value.Y.IsZero) or (Value.Z.IsZero) then Result.Infinity := true
else begin
     Result.Infinity := false;
     Result.SetCurveParams(Value.CurveParams, ComputeLambda);
     _Sqr_Fp6(Value.X, a);
     _Sqr_Fp6(Value.Y, b);
     _Sqr_Fp6(Value.Z, c);
     _Mul_Fp6(c, Value.CurveParams.BtwFp6, tmp);
     _Add_Fp6(tmp, tmp, D);
     _Add_Fp6(tmp, D, D);
     _Add_Fp6(Value.X, Value.Y, tmp);
     _Sqr_Fp6(tmp, E);
     _Sub_Fp6(E, a, E);
     _Sub_Fp6(E, b, E);
     _Add_Fp6(Value.Y, Value.Z, tmp);
     _Sqr_Fp6(tmp, F);
     _Sub_Fp6(F, b, F);
     _Sub_Fp6(F, c, F);
     _Add_Fp6(D, D, G);
     _Add_Fp6(D, G, G);
     _Sub_Fp6(b, G, tmp);
     _Mul_Fp6(tmp, E, Result.X);
     _Add_Fp6(b, G, tmp);
     _Sqr_Fp6(tmp, Result.Y);
     _Sqr_Fp6(D, tmp);
     _Add_Fp6(tmp, tmp, tmp);
     _Add_Fp6(tmp, tmp, tmp);
     _Sub_Fp6(Result.Y, tmp, Result.Y);
     _Add_Fp6(tmp, tmp, tmp);
     _Sub_Fp6(Result.Y, tmp, Result.Y);
     _Mul_Fp6(b, F, Result.Z);
     _Add_Fp6(Result.Z, Result.Z, Result.Z);
     _Add_Fp6(Result.Z, Result.Z, Result.Z);
     if Result.ComputeLigneValue then begin
                                      if Value.CurveParams.TwistMode=twDType then begin
                                      _Mul_Fp_Fp6(-Value.CurrentYp,F,Result.LineAtP36.a.a);
                                      _Sub_Fp6(d,b,Result.LineAtP36.a.b);
                                      _Add_Fp6(A,A,tmp);
                                      _Add_Fp6(A,tmp,tmp);
                                      _Mul_FP_Fp6(Value.CurrentXp,tmp,Result.LineAtP36.b.a);
                                      Result.LineAtP36.b.b.SetToZero;
                                      Result.LineAtP36.c.a.SetToZero;
                                      Result.LineAtP36.c.b.SetToZero;
                                      end
                                      else begin
                                           _Mul_Fp_Fp6(-Value.CurrentYp,F,Result.LineAtP36.a.b);
                                           _Sub_Fp6(d,b,Result.LineAtP36.a.a);
                                           _Add_Fp6(A,A,tmp);
                                           _Add_Fp6(A,tmp,tmp);
                                           _Mul_FP_Fp6(Value.CurrentXp,tmp,Result.LineAtP36.c.a);
                                           Result.LineAtP36.b.b.SetToZero;
                                           Result.LineAtP36.b.a.SetToZero;
                                           Result.LineAtP36.c.b.SetToZero;
                                           end;
                                      end;
  end;
end;

      { **********   Add two Fp6 Points (Projective coordinates)***************** }
      // Mixed addition (Left.Z is always equal to 1)
      // Costello et al's   https://cryptojedi.org/papers/edate-20100614.pdf
procedure _Add_Projective_Fp6_Point(Left {Q} , Right{T} : Fp6Point;var Result: Fp6Point; ComputeLambda: boolean = false);
var a, b, c, D,E,F,G,H,I,II,J,K,tmp1, tmp2,valx: Fp6Int;
begin
if Left.Infinity then Result := Right
else if Right.Infinity then Result := Left
else  begin
      if _Equals_Fp6(Left.X, Right.X) then begin
                                           if _Equals_Fp6(Left.Y, Right.Y) then _Double_Projective_Fp6_Point(Result, Result)
                                           else begin
                                                Result.Infinity := true;
                                                Result.Z.a.a := 1;
                                                Result.Z.a.b := 0;
                                                Result.Z.b.a := 0;
                                                Result.Z.b.b := 0;
                                                Result.SetCurveParams(Right.CurveParams, ComputeLambda);
                                                end;
                                           end
      else begin
          Result.SetCurveParams(Left.CurveParams, ComputeLambda);
          _Mul_Fp6(Right.Z, Left.X, a);
          _Sub_Fp6(Right.X, a, a);
          _Mul_Fp6(Right.Z, Left.Y, b);
          _Sub_Fp6(Right.Y, b, b);
          _Sqr_Fp6(a, tmp1);
          _Mul_Fp6(Right.X, tmp1, tmp2); // tmp2 used later as Right.X
          _Mul_Fp6(tmp1, a, c);
          _Mul_Fp6(b, b, tmp1);
          _Mul_Fp6(tmp1, Right.Z, D);
          _Add_Fp6(D, c, D);
          _Add_Fp6(tmp2, tmp2, tmp1);
          _Sub_Fp6(D, tmp1, D);
          _Sub_Fp6(tmp2, D, tmp2);
          _Mul_Fp6(b, tmp2, Result.Y);
          _Mul_Fp6(Right.Y, c, tmp2);
          _Sub_Fp6(Result.Y, tmp2, Result.Y);
          _Mul_Fp6(a, D, Result.X);
          _Mul_Fp6(Right.Z, c, Result.Z);
          if Result.ComputeLigneValue then begin
          if Left.CurveParams.TwistMode=twDType then begin
                                                           _Mul_FP_Fp6(-Right.CurrentYp,A,Result.LineAtP36.a.a);   ////L0,1
                                                           _Mul_Fp6(B,Left.X,tmp1);
                                                           _Mul_Fp6(A,Left.Y,tmp2);
                                                           _Sub_Fp6(tmp2,tmp1,Result.LineAtP36.a.b);
                                                           _Mul_FP_Fp6(Right.CurrentXp,B,Result.LineAtP36.b.a);
                                                           Result.LineAtP36.b.b.SetToZero;
                                                           Result.LineAtP36.c.a.SetToZero;
                                                           Result.LineAtP36.c.b.SetToZero;
                                                           end
          else begin
               _Mul_FP_Fp6(-Right.CurrentYp,A,Result.LineAtP36.a.b);   ////L0,1
               _Mul_Fp6(B,Left.X,tmp1);
               _Mul_Fp6(A,Left.Y,tmp2);
               _Sub_Fp6(tmp2,tmp1,Result.LineAtP36.a.a);
               _Mul_FP_Fp6(Right.CurrentXp,B,Result.LineAtP36.c.a);
               Result.LineAtP36.b.b.SetToZero;
               Result.LineAtP36.b.a.SetToZero;
               Result.LineAtP36.c.b.SetToZero;
               end;
               end;
      end
  end;
end;

      { **********   Add two Fp6 Points ***************** }
procedure _Add_Fp6_Point(Left, Right: Fp6Point; var Result: Fp6Point;CoordSys: CoordinatesSystem = csAffine);
begin
case CoordSys of
    csAffine:_Add_Affine_Fp6_Point(Left, Right, Result);
    csProjective:_Add_Projective_Fp6_Point(Left, Right, Result);
end;
end;

      { **********   Double an Fp6 Points ***************** }
procedure _Double_Fp6_Point(Value: Fp6Point; var Result: Fp6Point;
  CoordSys: CoordinatesSystem = csAffine);
begin
case CoordSys of
  csAffine:_Double_Affine_Fp6_Point(Value, Result);
  csProjective:_Double_Projective_Fp6_Point(Value, Result);
end;
end;

      { **********   Substract two Fp6 Points ***************** }
procedure _Sub_Fp6_Point(Left, Right: Fp6Point; var Result: Fp6Point;CoordSys: CoordinatesSystem = csAffine);
var  tmp: Fp6Point;
begin
_Neg_Fp6_Point(Right, tmp);
case CoordSys of
    csAffine:_Add_Affine_Fp6_Point(Left, tmp, Result);
    csProjective:_Add_Projective_Fp6_Point(Left, tmp, Result);
end;
end;

      { **********   Convert from Projective to Affine for an Fp6 Point ***** }
      // (X/Z,Y/Z)
procedure _Projective_To_Affine_Fp6Point(Value: Fp6Point; var Result: Fp6Point);
var  t1: Fp6Int;
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
      { **********   Convert from Projective to Affine for an Fp6 Point ***** }
      // (X/Z,Y/Z^2)
procedure _Projective2_To_Affine_Fp6Point(Value: Fp6Point; var Result: Fp6Point);
var  t1, t2: Fp6Int;
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


    { **********   Negation of an Fp6 Point ************* }
procedure _Neg_Fp6_Point(Value: Fp6Point; var Result: Fp6Point);
begin
Result := Value;
_Neg_Fp6(Result.Y, Result.Y);
end;


    { **********   Multiply a scalar with an Fp6 Point ******* }
procedure _Mul_Fp6_FpPoint(const Left: LInt; const Right: Fp6Point;var Result: Fp6Point; CoordSys: CoordinatesSystem = csAffine);
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
                         csAffine:_Double_Affine_Fp6_Point(Right, Result);
                         csProjective:_Double_Projective_Fp6_Point(Right, Result);
                        end;
                        end
       else for i := Left.BitLength - 1 downto 0 do begin
                                                    case CoordSys of
                                                    csAffine:_Double_Affine_Fp6_Point(Result,Result);
                                                    csProjective:_Double_Projective_Fp6_Point(Result, Result);
                                                    end;
                                                    if _Is_BitSet_At(Left,i) then begin
                                                                                  case CoordSys of
                                                                                    csAffine:_Add_Affine_Fp6_Point(Right, Result, Result);
                                                                                    csProjective:_Add_Projective_Fp6_Point(Right, Result, Result);
                                                                                  end;
                                                                                  end;
                                                    end;
       if CoordSys=csProjective then _Projective_To_Affine_Fp6Point(Result, Result);
       if _IsNeg(Left) then _Neg_Fp6_Point(Result,Result);
       end;
end;


      { **********   Multiply a scalar with an Fp6 Point (Negative Representation of the Exponent)******* }
procedure _Mul_NAF_Fp6_FpPoint(const Left: LInt; const Right: Fp6Point;var Result: Fp6Point; CoordSys: CoordinatesSystem = csAffine);
var  i: integer;
     loop: LIntArrayForm;
     NegLeft: Fp6Point;
begin
Result.CurveParams := Right.CurveParams;
_Neg_Fp6_Point(Right, NegLeft);
if (Right.Infinity) or (_IsNull(Left)) then Result.Infinity := true
else  begin
      Result.Infinity := true;
      if Left = 2 then _Double_Fp6_Point(Result, Result, CoordSys)
      else begin
           loop := Left.ToNafArray;
           for i := Length(loop) - 1 downto 0 do begin
                                                 _Double_Fp6_Point(Result, Result, CoordSys);
                                                 if loop[i] = 1 then _Add_Fp6_Point(Right, Result, Result, CoordSys);
                                                 if loop[i] = -1 then _Add_Fp6_Point(NegLeft, Result, Result, CoordSys);
                                                 end;
           _Projective_To_Affine_Fp6Point(Result, Result);
           end;
      end;
end;

      { **********   Compute Frobenius Map for an Fp6 Point ******* }
procedure _Frobenius_Map_i(const Value: Fp6Point; power:integer;var Result: Fp6Point);
var i:integer;
begin
if Value.Infinity then Result.Infinity:=true
else begin
     Result.SetCurveParams(Value.CurveParams);
     Result.X:=Value.X;
     Result.Y:=Value.Y;
     Result.Z:=Value.Z;
     for i:=0 to power-1 do begin
                          _Mul_Fp6(Result.CurveParams.FrobeniusMapConstX_Fp6,Result.X.toPowerP,Result.X);
                          _Mul_Fp6(Result.CurveParams.FrobeniusMapConstY_Fp6,Result.Y.toPowerP,Result.Y);
                          end;
    end;
end;
      { **********   Test if an Fp Point is on the Curve *********** }
function _Is_OnCurve_Fp6Point(Value: Fp6Point): boolean;
var tmp, tmp1: Fp6Int;
begin
if Value.Infinity then Result := true
else begin
     _Sqr_Fp6(Value.X, tmp);
     _Add_Fp6(tmp, Value.CurveParams.AtwFp6, tmp);
     _Mul_Fp6(tmp, Value.X, tmp1);
     _Add_Fp6(tmp1, Value.CurveParams.BtwFp6, tmp);
     _Sqr_Fp6(Value.Y, tmp1);
     Result:=(tmp.IsASquare) and (tmp1=tmp);
     end;
end;

      { **********   Test if two Fp6 Points are equal ******* }
function _Are_Equals_Fp6Points(Left, Right: Fp6Point): boolean;
begin
Result := (_Equals_Fp6(Left.X, Right.X)) and (_Equals_Fp6(Left.Y, Right.Y));
end;

      { **********   Evaluate an Fp6 Point from X value ****** }
procedure _Evaluate_Fp6_Point(Value: Fp6Int; var Result: Fp6Point);
var  tmp, tmp1: Fp6Int;
begin
Result.X := Value;
Result.Infinity := false;
_Sqr_Fp6(Result.X, tmp);
_Add_Fp6(tmp, Result.CurveParams.AtwFp6, tmp);
_Mul_Fp6(tmp, Result.X, tmp1);
_Add_Fp6(tmp1, Result.CurveParams.BtwFp6, tmp);
if tmp.IsASquare then  begin
                       _Sqrt_Fp6(tmp, Result.Y);
                       //_Neg_Fp6(Result.Y, tmp1);
                       //if _Equals_Fp6(Result.Y, tmp1) then Result.Y := tmp1;
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
procedure Fp6Point.SetCurveParams(Value: PtrCurveParams;
  ComputLmd: boolean = false);
begin
CurveParams := Value;
if CurveParams.TowerParam <> nil then LineAtP36.SetTowerParams(Value.TowerParam2);
X.SetTowerParams(Value.TowerParam);
Y.SetTowerParams(Value.TowerParam);
Z.SetTowerParams(Value.TowerParam);
ComputeLigneValue := ComputLmd;
Infinity := false;
end;

{ ******************************************************************************* }
{function Fp6Point.CompressToArray: TBytes;
var L:integer;
begin
L:=X.a.Field.p.Data.i32[-1]*4*4;
Setlength(Result,L);
Move(X.a.a.Data.i8[0],Result[0],Length(Result)shr 2);
Move(X.a.b.Data.i8[0],Result[L shr 2],Length(Result) shr 2);
Move(X.b.a.Data.i8[0],Result[L shr 1],Length(Result)shr 2);
Move(X.b.b.Data.i8[0],Result[L shr 1+L shr 2],Length(Result) shr 2);
end;}

{ ******************************************************************************* }
{procedure Fp6Point.DeCompressFromArray(a: TBytes);
begin
Move(a[0],X.a.a.Data.i8[0],Length(a) shr 2);
Move(a[Length(a) shr 2],X.a.b.Data.i8[0],Length(a) shr 2);
Move(a[Length(a) shr 1],X.b.a.Data.i8[0],Length(a) shr 2);
Move(a[Length(a) shr 1+length(a) shr 2],X.b.b.Data.i8[0],Length(a) shr 2);
X.a.a.Data.i16[-1]:=0;
X.a.b.Data.i16[-1]:=0;
X.b.a.Data.i16[-1]:=0;
X.b.b.Data.i16[-1]:=0;
X.a.a.Data.i16[-2]:=Length(a) shr 4;
X.a.b.Data.i16[-2]:=Length(a) shr 4;
X.b.a.Data.i16[-2]:=Length(a) shr 4;
X.b.b.Data.i16[-2]:=Length(a) shr 4;
SetPointValue(X);
end;}

{ ******************************************************************************* }
function Fp6Point.FrobeniusMap(pow:integer): Fp6Point;
begin
_Frobenius_Map_i(Self,pow, Result);
end;

{ ******************************************************************************* }
function Fp6Point.IsOnTheCurve: boolean;
begin
Result := _Is_OnCurve_Fp6Point(Self)
end;

{ ******************************************************************************* }
procedure Fp6Point.SetAsRandomPoint;
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
procedure Fp6Point.SetAsTorsionFromHash(hash: TBytes);
var i, j: integer;
    h1, h2: TBytes;
    tmp: Fp6Point;
begin
j:=0;
CurrentXp:=0;
CurrentYp:=0;
tmp.SetCurveParams (Self.CurveParams);
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
  X.b.a:=0;
  X.b.b:=0;
  X.c.a:=0;
  X.c.b:=0;
  SetPointValue(X);
  j:=j+1;
until not Infinity;
_Mul_Fp6_FpPoint(CurveParams.Htw,tmp,Self);
Self.Z.SetToZero;
Self.z.a.a:=1;
end;

{ ******************************************************************************* }
procedure Fp6Point.SetAsTorsionFromString(s:String);
begin
Self.SetAsTorsionFromHash(SHA256StringHash(s));
end;

{ ******************************************************************************* }
procedure Fp6Point.SetAsRandomTorsionPoint;
var tmp: Fp6Point;
begin
if CurveParams=nil then raise Exception.Create('Paramétre de courbe invalide');

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
                     //_Fast_Mul_Htw_Fp6Point(Self,tmp);
                     _Mul_Fp6_FpPoint(CurveParams.Htw,Self,tmp);
                     Self := tmp;
                     Self.Z.SetToZero;
                     Self.z.a.a:=1;
                     end;
until not Infinity;
{_Mul_Fp6_FpPoint(CurveParams.Rtw, Self, tmp, csaffine);   // For Testing Torsion point Validity
if not tmp.Infinity then showmessage('problem! G2');  }
end;

{ ******************************************************************************* }
procedure Fp6Point.SetPairingPointCoordinates(PtX, PtY: LInt);
begin
CurrentXp := PtX;
CurrentYp := PtY;
end;

{ ******************************************************************************* }
procedure Fp6Point.SetPointValue(Value: Fp6Int);
begin
_Evaluate_Fp6_Point(Value, Self);
end;

{ ******************************************************************************* }
procedure Fp6Point.SetToDefaultGenerator;
begin
Self.X:=CurveParams.KSS36TwistGeneratorX;
Self.Y:=CurveParams.KSS36TwistGeneratorY;
Self.Z.a.SetFormString('1');
Self.Z.b.SetFormString('0');
Self.Infinity:=false;
end;

{ ******************************************************************************* }
function Fp6Point.ToHexString: String;
begin
if Infinity then Result := 'Infinity'
else Result := '(' + X.ToHexString + ',' + Y.ToHexString + ')';
end;

{ ******************************************************************************* }
function Fp6Point.ToDecimalString: String;
begin
if Infinity then Result := 'Infinity'
else Result := '(' + X.ToDecimalString + ',' + Y.ToDecimalString + ')';
end;

{ ******************************************************************************* }
class operator Fp6Point.Add(Left, Right: Fp6Point): Fp6Point;
begin
_Add_Affine_Fp6_Point(Left, Right, Result);
end;

{ ******************************************************************************* }
class operator Fp6Point.Multiply(Left: LInt; Right { Q } : Fp6Point): Fp6Point;
begin
_Mul_Fp6_FpPoint(Left, Right, Result);
end;

{ ******************************************************************************* }
class operator Fp6Point.Subtract(Left, Right: Fp6Point): Fp6Point;
begin
_Sub_Fp6_Point(Left, Right, Result);
end;

{ ******************************************************************************* }
class operator Fp6Point.Negative(const Value: Fp6Point): Fp6Point;
begin
_Neg_Fp6_Point(Value, Result);
end;

{ ******************************************************************************* }
function Fp6Point.CompressToArray: TBytes;
var L:integer;
begin
if Infinity then begin
                 Setlength(result,0);
                 exit;
                 end;
L:=X.a.Field.p.Data.i32[-1]*4*6;
Setlength(Result,L);
Move(X.a.a.Data.i8[0],Result[0],Length(Result) div 6);
Move(X.a.b.Data.i8[0],Result[L div 6],Length(Result) div 6);
Move(X.b.a.Data.i8[0],Result[L div 6 *2],Length(Result) div 6);
Move(X.b.b.Data.i8[0],Result[L div 6 *3],Length(Result) div 6);
Move(X.c.a.Data.i8[0],Result[L div 6 *4],Length(Result) div 6);
Move(X.c.b.Data.i8[0],Result[L div 6 *5],Length(Result) div 6);
end;

{ ******************************************************************************* }
procedure Fp6Point.DeCompressFromArray(a: TBytes);
var L:integer;
begin
L:=length(a) div 6;
Move(a[0],X.a.a.Data.i8[0],L);
Move(a[L],X.a.b.Data.i8[0],L);
Move(a[L*2],X.b.a.Data.i8[0],L);
Move(a[L*3],X.b.b.Data.i8[0],L);
Move(a[L*4],X.c.a.Data.i8[0],L);
Move(a[L*5],X.c.b.Data.i8[0],L);
X.a.a.Data.i16[-1]:=0;
X.a.b.Data.i16[-1]:=0;
X.b.a.Data.i16[-1]:=0;
X.b.b.Data.i16[-1]:=0;
X.c.a.Data.i16[-1]:=0;
X.c.b.Data.i16[-1]:=0;
X.a.a.Data.i16[-2]:=L div 4;
X.a.b.Data.i16[-2]:=L div 4;
X.b.a.Data.i16[-2]:=L div 4;
X.b.b.Data.i16[-2]:=L div 4;
X.c.a.Data.i16[-2]:=L div 4;
X.c.b.Data.i16[-2]:=L div 4;
SetPointValue(X);
end;

class operator Fp6Point.Equal(const Left, Right: Fp6Point): boolean;
begin
Result := _Are_Equals_Fp6Points(Left, Right);
end;

{ ******************************************************************************* }
class operator Fp6Point.NotEqual(const Left, Right: Fp6Point): boolean;
begin
Result := not _Are_Equals_Fp6Points(Left, Right);
end;

end.

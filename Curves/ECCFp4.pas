unit ECCFp4;

interface

uses Fp4Arithmetic, System.SysUtils, LargeIntegers, Fp24Arithmetic,Fp16Arithmetic,Fp2Arithmetic,Fp8Arithmetic, HashFunctions, GeneralTypes,
     vcl.dialogs;
type
   //  {$M 204857600,2048586000}
  Fp4Point = record
    // Definition of a Point on EC Over Fp4

  public
    X, Y, Z: Fp4Int;
    ComputeLigneValue: boolean;
    CurveParams: PtrCurveParams;
    CurrentXp, CurrentYp: LInt;   // Used to store Coordinates of the E(Fp)'s point during Pairings
    LineAtP24: PFp24Int;             // Value of the line between the txo pairing points during Miller loop (On Fp24)
    LineAtP16: PFp16Int;             // Value of the line between the txo pairing points during Miller loop (On Fp16)
    Infinity: boolean;
    procedure SetCurveParams(Value: PtrCurveParams; ComputLmd: boolean = false);
    procedure SetPointValue(Value: Fp4Int);
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
    function FrobeniusMap(pow:integer): Fp4Point;
    class operator Add(Left, Right: Fp4Point): Fp4Point;
    class operator Subtract(Left, Right: Fp4Point): Fp4Point;
    class operator Multiply(Left: LInt; Right: Fp4Point): Fp4Point;
    class operator Equal(const Left, Right: Fp4Point): boolean;
    class operator NotEqual(const Left, Right: Fp4Point): boolean;
    class operator Negative(const Value: Fp4Point): Fp4Point;
  end;

  procedure _Add_Projective_Fp4_Point(Left { Q } , Right { T } : Fp4Point;  var Result: Fp4Point; ComputeLambda: boolean = false);
  procedure _Double_Projective_Fp4_Point(const Value: Fp4Point;  var Result: Fp4Point; ComputeLambda: boolean = false);
  procedure _Add_Affine_Fp4_Point(Left, Right: Fp4Point; var Result: Fp4Point;  ComputeLambda: boolean = false);
  procedure _Double_Affine_Fp4_Point(const Value: Fp4Point; var Result: Fp4Point;  ComputeLambda: boolean = false);
  procedure _Double_Fp4_Point(Value: Fp4Point; var Result: Fp4Point;  CoordSys: CoordinatesSystem = csAffine);
  procedure _Add_Fp4_Point(Left, Right: Fp4Point; var Result: Fp4Point;  CoordSys: CoordinatesSystem = csAffine);
  procedure _Sub_Fp4_Point(Left, Right: Fp4Point; var Result: Fp4Point;  CoordSys: CoordinatesSystem = csAffine);
  procedure _Evaluate_Fp4_Point(Value: Fp4Int; var Result: Fp4Point);
  procedure _Neg_Fp4_Point(Value: Fp4Point; var Result: Fp4Point);
  function _Is_OnCurve_Fp4Point(Value: Fp4Point): boolean;
  procedure _Frobenius_Map_i(const Value: Fp4Point; power:integer;var Result: Fp4Point);
  procedure _Mul_Fp4_FpPoint(const Left: LInt; const Right: Fp4Point;  var Result: Fp4Point; CoordSys: CoordinatesSystem = csAffine);
  procedure _Fast_Mul_Htw_Fp4Point(const Qx0: Fp4Point;var Result: Fp4Point);
  procedure _Projective_To_Affine_Fp4Point(Value: Fp4Point; var Result: Fp4Point);  //(X/Z,Y/Z)
  procedure _Projective2_To_Affine_Fp4Point(Value: Fp4Point; var Result: Fp4Point); //(X/Z,Y/Z^2)
  function _Are_Equals_Fp4Points(Left, Right: Fp4Point): boolean;


implementation

{ ******************************************************************************* }
            /// Procedures for Elliptic Curves Arithmetic over Fp4
{ ********************************************************************-*********** }

    { **********   Add two Fp4 Points  using Affine coordinates ***************** }
procedure _Add_Affine_Fp4_Point(Left, Right: Fp4Point; var Result: Fp4Point;ComputeLambda: boolean = false);
var  t: array [0 .. 4] of Fp4Int;
     Valx,Valy:Fp4Int;
begin
if Left.Infinity then Result := Right
else if Right.Infinity then Result := Left
else begin
     if _Equals_FP4(Left.X, Right.X) then begin
                                          if _Equals_FP4(Left.Y, Right.Y) then _Double_Affine_Fp4_Point(Left, Result)
                                          else begin
                                               Result.Z.a.a := 1;
                                               Result.Z.a.b := 0;
                                               Result.Z.b.a := 0;
                                               Result.Z.b.b := 0;
                                               Result.SetCurveParams(Right.CurveParams, ComputeLambda);
                                               Result.Infinity := true;
                                               end;
                                          end
     else begin
          Result.SetCurveParams(Left.CurveParams, ComputeLambda);
          ValX:=Right.X;
          ValY:=Right.Y;
          _Sub_FP4(Right.Y, Left.Y, t[0]);
          _Sub_FP4(Right.X, Left.X, t[1]);
          _Inv_FP4(t[1], t[2]);
          _Mul_FP4(t[0], t[2], t[4]); // Lambda
          _Sqr_FP4(t[4], t[3]);
          _Sub_FP4(t[3], Left.X, Result.X);
          _Sub_FP4(Result.X, Right.X, Result.X);
          _Sub_FP4(Left.X, Result.X, t[2]);
          _Mul_FP4(t[4], t[2], Result.Y);
          _Sub_FP4(Result.Y, Left.Y, Result.Y);
          if Result.ComputeLigneValue then begin
          if Left.CurveParams.TwistMode=twDType then begin
                                                     if Left.CurveParams.Family=cfKSS16 then begin
                                                                                             Result.LineAtP16.a.a.SetToZero;
                                                                                             _Sub_LInt(Right.CurveParams.FieldParam.p,Right.CurrentYp,Result.LineAtP16.a.a.a.a);
                                                                                             _Mul_FP4(t[4],Valx,t[3]);
                                                                                             _Sub_FP4(Valy,t[3],Result.LineAtP16.b.b);
                                                                                             _Mul_FP_FP4(Right.CurrentXp,t[4],Result.LineAtP16.b.a);
                                                                                             Result.LineAtP16.a.b.SetToZero;
                                                                                             end
                                                     else begin
                                                          Result.LineAtP24.a.a.SetToZero;
                                                          _Sub_LInt(Right.CurveParams.FieldParam.p,Right.CurrentYp,Result.LineAtP24.a.a.a.a);
                                                          _Mul_FP4(t[4],Valx,t[3]);
                                                          _Sub_FP4(Valy,t[3],Result.LineAtP24.a.b);
                                                          _Mul_FP_FP4(Right.CurrentXp,t[4],Result.LineAtP24.b.a);
                                                          Result.LineAtP24.b.b.SetToZero;
                                                          Result.LineAtP24.c.a.SetToZero;
                                                          Result.LineAtP24.c.b.SetToZero;
                                                          end
                                                    end
          else begin
           if Left.CurveParams.Family=cfKSS16 then begin
                                                    Result.LineAtP16.b.b.SetToZero;
                                                    _Sub_LInt(Right.CurveParams.FieldParam.p,Right.CurrentYp,Result.LineAtP16.b.b.a.a);
                                                    _Mul_FP4(t[4],Valx,t[3]);
                                                    _Sub_FP4(Valy,t[3],Result.LineAtP16.a.a);
                                                    _Mul_FP_FP4(Right.CurrentXp,t[4],Result.LineAtP16.a.b);
                                                    Result.LineAtP16.b.a.SetToZero;
                                                    end
           else begin
                 Result.LineAtP24.a.b.SetToZero;
                 _Sub_LInt(Right.CurveParams.FieldParam.p,Right.CurrentYp,Result.LineAtP24.a.b.a.a);
                 _Mul_FP4(t[4],Valx,t[3]);
                 _Sub_FP4(Valy,t[3],Result.LineAtP24.a.a);
                 _Mul_FP_FP4(Right.CurrentXp,t[4],Result.LineAtP24.c.a);
                 Result.LineAtP24.b.b.SetToZero;
                 Result.LineAtP24.b.a.SetToZero;
                 Result.LineAtP24.c.b.SetToZero;
                 end;
              end;
          end;
          end;
     end;
end;

      { **********   Double an Fp4 Point using Affine coordinates ***************** }
procedure _Double_Affine_Fp4_Point(const Value: Fp4Point; var Result: Fp4Point;ComputeLambda: boolean = false);
var t: array [0 .. 3] of Fp4Int;
    ValX, ValY: Fp4Int;
begin
if (Value.Infinity) or (Value.Y.IsZero) then Result.Infinity := true
else  begin
      Result.Infinity := false;
      Result.SetCurveParams(Value.CurveParams, ComputeLambda);
      ValX := Value.X;
      ValY := Value.Y;
      _Sqr_FP4(Value.X, t[0]);
      _Mul_FP_FP4(LInt(3), t[0], t[1]);
       _Add_FP4(t[1],Value.CurveParams.AtwFp4 ,t[1]);   /// Atw is Zero on BLS curves
      _Mul_FP_FP4(LInt(2), Value.Y, t[3]);
      _Inv_FP4(t[3], t[3]);
      _Mul_FP4(t[1], t[3], t[2]); // Lambda
      _Sqr_FP4(t[2], t[3]);
      _Sub_FP4(t[3], Value.X, t[3]);
      _Sub_FP4(t[3], Value.X, Result.X);
      _Sub_FP4(ValX, Result.X, t[3]);
      _Mul_FP4(t[2], t[3], Result.Y);
      _Sub_FP4(Result.Y, ValY, Result.Y);
      if Result.ComputeLigneValue then begin
      if Value.CurveParams.TwistMode=twDType then begin
                                                  if Value.CurveParams.Family=cfKSS16 then begin
                                                                                           Result.LineAtP16.a.a.SetToZero;
                                                                                           _Sub_LInt(Value.CurveParams.FieldParam.p,Value.CurrentYp,Result.LineAtP16.a.a.a.a);
                                                                                           _Mul_FP4(t[2],Valx,t[3]);
                                                                                           _Sub_FP4(Valy,t[3],Result.LineAtP16.b.b);
                                                                                           _Mul_FP_FP4(Value.CurrentXp,t[2],Result.LineAtP16.b.a);
                                                                                           Result.LineAtP16.a.b.SetToZero;
                                                                                           end
                                                  else begin
                                                       Result.LineAtP24.a.a.SetToZero;
                                                       _Sub_LInt(Value.CurveParams.FieldParam.p,Value.CurrentYp,Result.LineAtP24.a.a.a.a);
                                                       _Mul_FP4(t[2],Valx,t[3]);
                                                       _Sub_FP4(Valy,t[3],Result.LineAtP24.a.b);
                                                       _Mul_FP_FP4(Value.CurrentXp,t[2],Result.LineAtP24.b.a);
                                                       Result.LineAtP24.b.b.SetToZero;
                                                       Result.LineAtP24.c.a.SetToZero;
                                                       Result.LineAtP24.c.b.SetToZero;
                                                       end;
                                                  end
      else begin
           if Value.CurveParams.Family=cfKSS16 then begin
                                                    Result.LineAtP16.b.b.SetToZero;
                                                    _Sub_LInt(Value.CurveParams.FieldParam.p,Value.CurrentYp,Result.LineAtP16.b.b.a.a);
                                                    _Mul_FP4(t[2],Valx,t[3]);
                                                    _Sub_FP4(Valy,t[3],Result.LineAtP16.a.a);
                                                    _Mul_FP_FP4(Value.CurrentXp,t[2],Result.LineAtP16.a.b);
                                                    Result.LineAtP16.b.a.SetToZero;
                                                    end
           else begin
                Result.LineAtP24.a.b.SetToZero;
                _Sub_LInt(Value.CurveParams.FieldParam.p,Value.CurrentYp,Result.LineAtP24.a.b.a.a);
                _Mul_FP4(t[2],Valx,t[3]);
                _Sub_FP4(Valy,t[3],Result.LineAtP24.a.a);
                _Mul_FP_FP4(Value.CurrentXp,t[2],Result.LineAtP24.c.a);
                Result.LineAtP24.b.b.SetToZero;
                Result.LineAtP24.b.a.SetToZero;
                Result.LineAtP24.c.b.SetToZero;
                end;
           end;
      end;
      end;
end;

      { **********   Double an Fp4 Point (Projective coordinates (X/Z,Y/Z))***************** }
      // Costello et al's   https://cryptojedi.org/papers/edate-20100614.pdf
procedure _Double_Projective_Fp4_Point(const Value: Fp4Point;var Result: Fp4Point; ComputeLambda: boolean = false);
var tmp,L00,L01,L10, a, b, c, D, E, F, G: Fp4Int;
begin
if (Value.Infinity) or (Value.Y.IsZero) or (Value.Z.IsZero) then Result.Infinity := true
else begin
     Result.Infinity := false;
     Result.SetCurveParams(Value.CurveParams, ComputeLambda);
     if  Value.CurveParams.AtwFp4.IsZero then begin    //// This part is for Sexstic Twists of BLS24 curves with equation y^2=x^3+b
                                             _Sqr_FP4(Value.X, a);
                                             _Sqr_FP4(Value.Y, b);
                                             _Sqr_FP4(Value.Z, c);
                                             _Mul_FP4(c, Value.CurveParams.BtwFp4, tmp);
                                             _Add_FP4(tmp, tmp, D);
                                             _Add_FP4(tmp, D, D);
                                             _Add_FP4(Value.X, Value.Y, tmp);
                                             _Sqr_FP4(tmp, E);
                                             _Sub_FP4(E, a, E);
                                             _Sub_FP4(E, b, E);
                                             _Add_FP4(Value.Y, Value.Z, tmp);
                                             _Sqr_FP4(tmp, F);
                                             _Sub_FP4(F, b, F);
                                             _Sub_FP4(F, c, F);
                                             _Add_FP4(D, D, G);
                                             _Add_FP4(D, G, G);
                                             _Sub_FP4(b, G, tmp);
                                             _Mul_FP4(tmp, E, Result.X);
                                             _Add_FP4(b, G, tmp);
                                             _Sqr_FP4(tmp, Result.Y);
                                             _Sqr_FP4(D, tmp);
                                             _Add_FP4(tmp, tmp, tmp);
                                             _Add_FP4(tmp, tmp, tmp);
                                             _Sub_FP4(Result.Y, tmp, Result.Y);
                                             _Add_FP4(tmp, tmp, tmp);
                                             _Sub_FP4(Result.Y, tmp, Result.Y);
                                             _Mul_FP4(b, F, Result.Z);
                                             _Add_FP4(Result.Z, Result.Z, Result.Z);
                                             _Add_FP4(Result.Z, Result.Z, Result.Z);
                                             if Result.ComputeLigneValue then begin
                                                                               if Value.CurveParams.TwistMode=twDType then begin
                                                                                                                             _Mul_Fp_FP4(-Value.CurrentYp,F,Result.LineAtP24.a.a);
                                                                                                                             _Sub_FP4(d,b,Result.LineAtP24.a.b);
                                                                                                                             _Nomod_Add_FP4(A,A,tmp);
                                                                                                                             _Nomod_Add_FP4(A,tmp,tmp);
                                                                                                                             _Mul_FP_FP4(Value.CurrentXp,tmp,Result.LineAtP24.b.a);
                                                                                                                             Result.LineAtP24.b.b.SetToZero;
                                                                                                                             Result.LineAtP24.c.a.SetToZero;
                                                                                                                             Result.LineAtP24.c.b.SetToZero;
                                                                                                                             end
                                                                               else begin
                                                                                _Mul_Fp_FP4(-Value.CurrentYp,F,Result.LineAtP24.a.b);
                                                                                _Sub_FP4(d,b,Result.LineAtP24.a.a);
                                                                                _Nomod_Add_FP4(A,A,tmp);
                                                                                _Nomod_Add_FP4(A,tmp,tmp);
                                                                                _Mul_FP_FP4(Value.CurrentXp,tmp,Result.LineAtP24.c.a);
                                                                                Result.LineAtP24.b.b.SetToZero;
                                                                                Result.LineAtP24.b.a.SetToZero;
                                                                                Result.LineAtP24.c.b.SetToZero;
                                                                                end;
                                                                              end;
                                            end
      else begin         //// This part is for Quartic Twists of KSS16 curves with equation y^2=x^3+a*x
           _Sqr_FP4(Value.X, a);
           _Sqr_FP4(Value.Y, b);
           _Sqr_FP4(Value.Z, c);
           _Mul_FP4(c, Value.CurveParams.AtwFp4, D);
           _Sub_FP4(a, d, G);    //G is used as temporary variable only
           _Add_FP4(A,A,L10);
           _Add_FP4(A,L10,L10);
           _Add_FP4(L10,d,L10);
           _Add_FP4(L10,L10,L10);
           _Mul_FP4(Value.Z,L10,L10);
           _Add_FP4(Value.Y,Value.Z,L01);
           _Sqr_FP4(L01,L01);
           _Sub_FP4(L01,B,L01);
           _Sub_FP4(L01,C,L01);
           _Add_FP4(L01,L01,L01);
           _Sub_FP4(A,D,L00);
           _Add_FP4(L00,Value.X,L00);
           _Sqr_FP4(L00,L00);
           _Sqr_FP4(G,Result.X);
           _Sub_FP4(L00,Result.X,L00);
           _Sub_FP4(L00,A,L00);
           _Add_FP4(a, d, tmp);
           _Sqr_FP4(tmp,E);
           _Add_FP4(E,E,E);
           _Sub_FP4(E,Result.X,E);
           _Add_FP4(G,Value.Y,G);
           _Sqr_FP4(G,G);
           _Sub_FP4(G,B,G);
           _Sub_FP4(G,Result.X,F);
           _Mul_FP4(E,F,Result.Y);
           _Add_FP4(B,B,Result.Z);
           _Add_FP4(Result.Z,Result.Z,Result.Z);
           if Result.ComputeLigneValue then begin
                                            if Value.CurveParams.TwistMode=twDType then begin
                                                                                             _Mul_Fp_FP4(Value.CurrentYp,L01,Result.LineAtP16.a.a);
                                                                                             Result.LineAtP16.b.b:=L00;
                                                                                             _Mul_FP_FP4(-Value.CurrentXp,L10,Result.LineAtP16.b.a);
                                                                                             Result.LineAtP16.a.b.SetToZero;
                                                                                             end
                                            else begin
                                                  _Mul_Fp_FP4(Value.CurrentYp,L01,Result.LineAtP16.b.b);
                                                  Result.LineAtP16.a.a:=L00;
                                                  _Mul_FP_FP4(-Value.CurrentXp,L10,Result.LineAtP16.a.b);
                                                  Result.LineAtP16.b.a.SetToZero;
                                                 end;
                                             end;
       end;
  end;
end;

      { **********   Add two Fp4 Points (Projective coordinates)***************** }
      // Mixed addition (Left.Z is always equal to 1)
      // Costello et al's   https://cryptojedi.org/papers/edate-20100614.pdf
procedure _Add_Projective_Fp4_Point(Left {Q} , Right{T} : Fp4Point;var Result: Fp4Point; ComputeLambda: boolean = false);
var a, b, c, D,E,F,G,H,I,II,J,K,tmp1, tmp2,valx: Fp4Int;
begin
if Left.Infinity then Result := Right
else if Right.Infinity then Result := Left
else  begin
      if _Equals_FP4(Left.X, Right.X) then begin
                                           if _Equals_FP4(Left.Y, Right.Y) then _Double_Projective_Fp4_Point(Result, Result)
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
     if  Left.CurveParams.AtwFp4.IsZero then begin    //// This part is for Sexstic Twists of BLS24 curves with equation y^2=x^3+b
                                             _Mul_FP4(Right.Z, Left.X, a);
                                             _Sub_FP4(Right.X, a, a);
                                             _Mul_FP4(Right.Z, Left.Y, b);
                                             _Sub_FP4(Right.Y, b, b);
                                             _Sqr_FP4(a, tmp1);
                                             _Mul_FP4(Right.X, tmp1, tmp2); // tmp2 used later as Right.X
                                             _Mul_FP4(tmp1, a, c);
                                             _Mul_FP4(b, b, tmp1);
                                             _Mul_FP4(tmp1, Right.Z, D);
                                             _Add_FP4(D, c, D);
                                             _Add_FP4(tmp2, tmp2, tmp1);
                                             _Sub_FP4(D, tmp1, D);
                                             _Sub_FP4(tmp2, D, tmp2);
                                             _Mul_FP4(b, tmp2, Result.Y);
                                             _Mul_FP4(Right.Y, c, tmp2);
                                             _Sub_FP4(Result.Y, tmp2, Result.Y);
                                             _Mul_FP4(a, D, Result.X);
                                             _Mul_FP4(Right.Z, c, Result.Z);
                                             if Result.ComputeLigneValue then begin
                                                                          if Left.CurveParams.TwistMode=twDType then
                                                                                          begin
                                                                                          _Mul_FP_FP4(-Right.CurrentYp,A,Result.LineAtP24.a.a);   ////L0,1
                                                                                          _Mul_FP4(B,Left.X,tmp1);
                                                                                          _Mul_FP4(A,Left.Y,tmp2);
                                                                                          _Sub_FP4(tmp2,tmp1,Result.LineAtP24.a.b);
                                                                                          _Mul_FP_FP4(Right.CurrentXp,B,Result.LineAtP24.b.a);
                                                                                          Result.LineAtP24.b.b.SetToZero;
                                                                                          Result.LineAtP24.c.a.SetToZero;
                                                                                          Result.LineAtP24.c.b.SetToZero;
                                                                                          end
                                                                         else begin
                                                                              _Mul_FP_FP4(-Right.CurrentYp,A,Result.LineAtP24.a.b);   ////L0,1
                                                                              _Mul_FP4(B,Left.X,tmp1);
                                                                              _Mul_FP4(A,Left.Y,tmp2);
                                                                              _Sub_FP4(tmp2,tmp1,Result.LineAtP24.a.a);
                                                                              _Mul_FP_FP4(Right.CurrentXp,B,Result.LineAtP24.c.a);
                                                                              Result.LineAtP24.b.b.SetToZero;
                                                                              Result.LineAtP24.b.a.SetToZero;
                                                                              Result.LineAtP24.c.b.SetToZero;
                                                                              end;
                                                                            end;
                                           end
     else begin   //// This part is for Quartic Twists of KSS16 curves with equation y^2=x^3+a*x
            valx:=Right.X;                /// Can optimize Further since Left.Z=1 Always
           _Sqr_FP4(Right.Z,a);
           _Add_FP4(Right.Z,Right.Z,c);
           _Mul_FP4(Left.X,Right.Z,E);
           _Mul_FP4(Left.Y,A,G);
           _Sub_FP4(Right.X,E,H);
           _Sub_FP4(Right.Y,G,I);
           _Add_FP4(I,I,I);
           _Sqr_FP4(I,II);
           _Mul_FP4(Right.Z,H,J);
           _Add_FP4(J,J,J);
           _Mul_FP4(j,H,K);
           _Add_FP4(K,K,K);
           _Add_FP4(K,K,K);
           _Add_FP4(ValX,E,Result.X);
           _Mul_FP4(Result.X,K,Result.X);
           _Neg_FP4(Result.X,Result.X);
           _Add_FP4(Result.X,II,Result.X);
           _Add_FP4(Result.X,II,Result.X);
           _Sqr_FP4(J,Result.Z);
           _Sqr_FP4(K,tmp2);
           _Mul_FP4(tmp2,Right.Y,tmp2);
           _Add_FP4(J,I,tmp1);
           _Sqr_FP4(tmp1,tmp1);
           _Sub_FP4(tmp1,Result.Z,tmp1);
           _Sub_FP4(tmp1,II,tmp1);
           _Mul_FP4(ValX,K,Result.Y);
           _Sub_FP4(Result.Y,Result.X,Result.Y);
           _Mul_FP4(Result.Y,tmp1,Result.Y);
           _Sub_FP4(Result.Y,tmp2,Result.Y);
           _Add_FP4(Result.Z,Result.Z,Result.Z);
          if Result.ComputeLigneValue then begin
                                            if Left.CurveParams.TwistMode=twDType then begin
                                                                                             _Mul_FP_FP4(Right.CurrentYp,J,Result.LineAtP16.a.a);   ////L0,1
                                                                                             _Mul_FP4(I,Left.X,tmp1);
                                                                                             _Mul_FP4(J,Left.Y,tmp2);
                                                                                             _Sub_FP4(tmp1,tmp2,Result.LineAtP16.b.b);
                                                                                             _Mul_FP_FP4(-Right.CurrentXp,I,Result.LineAtP16.b.a);
                                                                                             Result.LineAtP16.a.b.SetToZero;
                                                                                             end
                                            else begin
                                                 _Mul_FP_FP4(Right.CurrentYp,J,Result.LineAtP16.b.b);   ////L0,1
                                                 _Mul_FP4(I,Left.X,tmp1);
                                                 _Mul_FP4(J,Left.Y,tmp2);
                                                 _Sub_FP4(tmp1,tmp2,Result.LineAtP16.a.a);
                                                 _Mul_FP_FP4(-Right.CurrentXp,I,Result.LineAtP16.a.b);
                                                 Result.LineAtP16.b.a.SetToZero;
                                                 end;
                                         end;
     end;
   end;
  end;
end;

      { **********   Add two Fp4 Points ***************** }
procedure _Add_Fp4_Point(Left, Right: Fp4Point; var Result: Fp4Point;CoordSys: CoordinatesSystem = csAffine);
begin
case CoordSys of
    csAffine:_Add_Affine_Fp4_Point(Left, Right, Result);
    csProjective:_Add_Projective_Fp4_Point(Left, Right, Result);
end;
end;

      { **********   Double an Fp4 Points ***************** }
procedure _Double_Fp4_Point(Value: Fp4Point; var Result: Fp4Point;
  CoordSys: CoordinatesSystem = csAffine);
begin
case CoordSys of
  csAffine:_Double_Affine_Fp4_Point(Value, Result);
  csProjective:_Double_Projective_Fp4_Point(Value, Result);
end;
end;

      { **********   Substract two Fp4 Points ***************** }
procedure _Sub_Fp4_Point(Left, Right: Fp4Point; var Result: Fp4Point;CoordSys: CoordinatesSystem = csAffine);
var  tmp: Fp4Point;
begin
_Neg_Fp4_Point(Right, tmp);
case CoordSys of
    csAffine:_Add_Affine_Fp4_Point(Left, tmp, Result);
    csProjective:_Add_Projective_Fp4_Point(Left, tmp, Result);
end;
end;

      { **********   Convert from Projective to Affine for an Fp4 Point ***** }
      // (X/Z,Y/Z)
procedure _Projective_To_Affine_Fp4Point(Value: Fp4Point; var Result: Fp4Point);
var  t1: Fp4Int;
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
      { **********   Convert from Projective to Affine for an Fp4 Point ***** }
      // (X/Z,Y/Z^2)
procedure _Projective2_To_Affine_Fp4Point(Value: Fp4Point; var Result: Fp4Point);
var  t1, t2: Fp4Int;
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


    { **********   Negation of an Fp4 Point ************* }
procedure _Neg_Fp4_Point(Value: Fp4Point; var Result: Fp4Point);
begin
Result := Value;
_Neg_FP4(Result.Y, Result.Y);
end;


    { **********   Multiply a scalar with an Fp4 Point ******* }
procedure _Mul_Fp4_FpPoint(const Left: LInt; const Right: Fp4Point;var Result: Fp4Point; CoordSys: CoordinatesSystem = csAffine);
// Right Should be different than Result
var  i: integer;
begin
Result.CurveParams := Right.CurveParams;
if  Right.CurveParams<>nil then Result.CurveParams:=Right.CurveParams;
if (Right.Infinity) or (_IsNull(Left)) then Result.Infinity := true
  else begin
       Result.Infinity := true;
       if Left = 2 then begin
                        case CoordSys of
                         csAffine:_Double_Affine_Fp4_Point(Right, Result);
                         csProjective:_Double_Projective_Fp4_Point(Right, Result);
                        end;
                        end
       else for i := Left.BitLength - 1 downto 0 do begin
                                                    case CoordSys of
                                                    csAffine:_Double_Affine_Fp4_Point(Result,Result);
                                                    csProjective:_Double_Projective_Fp4_Point(Result, Result);
                                                    end;
                                                    if _Is_BitSet_At(Left,i) then begin
                                                                                  case CoordSys of
                                                                                    csAffine:_Add_Affine_Fp4_Point(Right, Result, Result);
                                                                                    csProjective:_Add_Projective_Fp4_Point(Right, Result, Result);
                                                                                  end;
                                                                                  end;
                                                    end;
       if CoordSys=csProjective then begin
                                     if Right.CurveParams.Family=cfKSS16 then _Projective2_To_Affine_Fp4Point(Result, Result) /// in this case the form (X/Z,Y/Z^2) is better !
                                     else _Projective_To_Affine_Fp4Point(Result, Result);
                                     end;
       if _IsNeg(Left) then _Neg_Fp4_Point(Result,Result);
       end;
end;

      { **********   Multiply the Cofactor with an Fp4 Point ******* }
      // Fast Multiplication of A point by the Twist Cofactor (Generated By  Luis Dominguez)
      // https://eprint.iacr.org/2008/530.pdf
procedure _Fast_Mul_Htw_Fp4Point(const Qx0: Fp4Point;var Result: Fp4Point);
var  xA,xB,xC,t0,t1,Qx0_,Qx1,Qx1_,Qx2,Qx2_,Qx3,Qx3_,
     Qx4,Qx4_,Qx5,Qx5_,Qx6,Qx6_,Qx7,Qx7_,Qx8,Qx8_: Fp4Point;
begin
Result.CurveParams:=Qx0.CurveParams;
if Qx0.CurveParams<>nil then  Result.CurveParams:=Qx0.CurveParams;
_Neg_Fp4_Point(Qx0,Qx0_);
_Mul_Fp4_FpPoint(Qx0.CurveParams.u,Qx0,Qx1,csProjective);
_Neg_Fp4_Point(Qx1,Qx1_);
_Mul_Fp4_FpPoint(Qx0.CurveParams.u,Qx1,Qx2,csProjective);
_Neg_Fp4_Point(Qx2,Qx2_);
_Mul_Fp4_FpPoint(Qx0.CurveParams.u,Qx2,Qx3,csProjective);
_Neg_Fp4_Point(Qx3,Qx3_);
_Mul_Fp4_FpPoint(Qx0.CurveParams.u,Qx3,Qx4,csProjective);
_Neg_Fp4_Point(Qx4,Qx4_);
_Mul_Fp4_FpPoint(Qx0.CurveParams.u,Qx4,Qx5,csProjective);
_Neg_Fp4_Point(Qx5,Qx5_);
_Mul_Fp4_FpPoint(Qx0.CurveParams.u,Qx5,Qx6,csProjective);
_Neg_Fp4_Point(Qx6,Qx6_);
_Mul_Fp4_FpPoint(Qx0.CurveParams.u,Qx6,Qx7,csProjective);
_Neg_Fp4_Point(Qx7,Qx7_);
_Mul_Fp4_FpPoint(Qx0.CurveParams.u,Qx7,Qx8,csProjective);
_Neg_Fp4_Point(Qx8,Qx8_);
xA:=Qx0;
xC:=Qx7;
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx2,4,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
xB:=Qx0;
xC:=Qx7;
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx2,4,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Add_Affine_Fp4_Point(xA,xB,t0); //
xB:=Qx2_;
xC:=Qx3_;
_Add_Affine_Fp4_Point(xB,xC,xB);
xC:=Qx8_;
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx2,1,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx3_,1,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx1,6,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Add_Affine_Fp4_Point(t0,xB,t0);     //
xB:=Qx4;
xC:=Qx5_;
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx0_,4,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx4_,4,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx0,5,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx1_,5,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx2_,5,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx3,5,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Add_Affine_Fp4_Point(t0,xB,t0); //
xA:=Qx1;
_Frobenius_Map_i(Qx0_,1,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx1,1,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx4_,1,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx5,1,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx0,2,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx1_,2,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx4_,2,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx5,2,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx0,3,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx1_,3,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx4_,3,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx5,3,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx1,4,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx3,4,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx0_,6,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
_Frobenius_Map_i(Qx2_,6,xC);
_Add_Affine_Fp4_Point(xA,xC,xA);
xB:=Qx4;
xC:=Qx5_;
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx0_,4,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx4_,4,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx0,5,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx1_,5,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx2_,5,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Frobenius_Map_i(Qx3,5,xC);
_Add_Affine_Fp4_Point(xB,xC,xB);
_Add_Affine_Fp4_Point(xA,xB,t1);
_Double_Affine_Fp4_Point(t0,t0);
_Add_Affine_Fp4_Point(t0,t1,Result);
end;

      { **********   Multiply a scalar with an Fp4 Point (Negative Representation of the Exponent)******* }
procedure _Mul_NAF_Fp4_FpPoint(const Left: LInt; const Right: Fp4Point;var Result: Fp4Point; CoordSys: CoordinatesSystem = csAffine);
var  i: integer;
     loop: LIntArrayForm;
     NegLeft: Fp4Point;
begin
Result.CurveParams := Right.CurveParams;
_Neg_Fp4_Point(Right, NegLeft);
if (Right.Infinity) or (_IsNull(Left)) then Result.Infinity := true
else  begin
      Result.Infinity := true;
      if Left = 2 then _Double_Fp4_Point(Result, Result, CoordSys)
      else begin
           loop := Left.ToNafArray;
           for i := Length(loop) - 1 downto 0 do begin
                                                 _Double_Fp4_Point(Result, Result, CoordSys);
                                                 if loop[i] = 1 then _Add_Fp4_Point(Right, Result, Result, CoordSys);
                                                 if loop[i] = -1 then _Add_Fp4_Point(NegLeft, Result, Result, CoordSys);
                                                 end;
           _Projective_To_Affine_Fp4Point(Result, Result);
           end;
      end;
end;

      { **********   Compute Frobenius Map for an Fp4 Point ******* }
procedure _Frobenius_Map_i(const Value: Fp4Point; power:integer;var Result: Fp4Point);
var i:integer;
begin
if Value.Infinity then Result.Infinity:=true
else begin
     Result.SetCurveParams(Value.CurveParams);
     Result.X:=Value.X;
     Result.Y:=Value.Y;
     Result.Z:=Value.Z;
     for i:=0 to power-1 do begin
                          _Mul_Fp4(Result.CurveParams.FrobeniusMapConstX_FP4,Result.X.toPowerP,Result.X);
                          _Mul_Fp4(Result.CurveParams.FrobeniusMapConstY_Fp4,Result.Y.toPowerP,Result.Y);
                          end;
    end;
end;
      { **********   Test if an Fp Point is on the Curve *********** }
function _Is_OnCurve_Fp4Point(Value: Fp4Point): boolean;
var tmp, tmp1: Fp4Int;
begin
if Value.Infinity then Result := true
else begin
     _Sqr_FP4(Value.X, tmp);
     _Add_FP4(tmp, Value.CurveParams.AtwFp4, tmp);
     _Mul_FP4(tmp, Value.X, tmp1);
     _Add_FP4(tmp1, Value.CurveParams.BtwFp4, tmp);
     _Sqr_FP4(Value.Y, tmp1);
     Result:=(tmp.IsASquare) and (tmp1=tmp);
     end;
end;

      { **********   Test if two Fp4 Points are equal ******* }
function _Are_Equals_Fp4Points(Left, Right: Fp4Point): boolean;
begin
Result := (_Equals_FP4(Left.X, Right.X)) and (_Equals_FP4(Left.Y, Right.Y));
end;

      { **********   Evaluate an Fp4 Point from X value ****** }
procedure _Evaluate_Fp4_Point(Value: Fp4Int; var Result: Fp4Point);
var  tmp, tmp1: Fp4Int;
begin
Result.X := Value;
Result.Infinity := false;
_Sqr_FP4(Result.X, tmp);
_Add_FP4(tmp, Result.CurveParams.AtwFp4, tmp);
_Mul_FP4(tmp, Result.X, tmp1);
_Add_FP4(tmp1, Result.CurveParams.BtwFp4, tmp);
if tmp.IsASquare then  begin
                       _Sqrt_FP4(tmp, Result.Y);
                       //_Neg_FP4(Result.Y, tmp1);
                       //if _Equals_FP4(Result.Y, tmp1) then Result.Y := tmp1;
                       Result.Z.a.a := 1;
                       Result.Z.a.b := 0;
                       Result.Z.b.a := 0;
                       Result.Z.b.b := 0;
                       end
else Result.Infinity := true;
end;

{ ******************************************************************************* }
          // Definitions of an Fp4 Point operators and functions
{ ******************************************************************************* }

procedure Fp4Point.SetCurveParams(Value: PtrCurveParams;
  ComputLmd: boolean = false);
begin
CurveParams := Value;
if CurveParams.TowerParam <> nil then begin
                                      if CurveParams.Family=cfBLS24 then begin
                                                                         new(LineAtP24);
                                                                         LineAtP24.SetTowerParams(Value.TowerParam2);
                                                                         end
                                      else if CurveParams.Family=cfKSS16 then begin
                                                                              new(LineAtP16);
                                                                              LineAtP16.SetTowerParams(Value.TowerParam2);
                                                                              end;
                                      end;
X.SetTowerParams(Value.TowerParam);
Y.SetTowerParams(Value.TowerParam);
Z.SetTowerParams(Value.TowerParam);
ComputeLigneValue := ComputLmd;
Infinity := false;
end;

{ ******************************************************************************* }
function Fp4Point.CompressToArray: TBytes;
var L:integer;
begin
L:=X.a.Field.p.Data.i32[-1]*4*4;
Setlength(Result,L);
Move(X.a.a.Data.i8[0],Result[0],Length(Result)shr 2);
Move(X.a.b.Data.i8[0],Result[L shr 2],Length(Result) shr 2);
Move(X.b.a.Data.i8[0],Result[L shr 1],Length(Result)shr 2);
Move(X.b.b.Data.i8[0],Result[L shr 1+L shr 2],Length(Result) shr 2);
end;

{ ******************************************************************************* }
procedure Fp4Point.DeCompressFromArray(a: TBytes);
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
end;

{ ******************************************************************************* }
function Fp4Point.FrobeniusMap(pow:integer): Fp4Point;
begin
_Frobenius_Map_i(Self,pow, Result);
end;

{ ******************************************************************************* }
function Fp4Point.IsOnTheCurve: boolean;
begin
Result := _Is_OnCurve_Fp4Point(Self)
end;

{ ******************************************************************************* }
procedure Fp4Point.SetAsRandomPoint;
begin
Infinity := false;
Randomize;
repeat
    X.a.SetToRandom;
    X.b.SetToRandom;
    SetPointValue(X);
until not Infinity;
end;

{ ******************************************************************************* }
procedure Fp4Point.SetAsTorsionFromHash(hash: TBytes);
var i, j: integer;
    h1, h2: TBytes;
    tmp: Fp4Point;
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
  X.b.a:=0;
  X.b.b:=0;
  SetPointValue(X);
  j:=j+1;
until not Infinity;
if CurveParams.Family=cfKSS16 then _Mul_Fp4_FpPoint(CurveParams.htw, tmp, Self, csaffine)
else _Fast_Mul_Htw_Fp4Point(tmp,Self);
Self.Z.SetToZero;
Self.z.a.a:=1;
end;

{ ******************************************************************************* }
procedure Fp4Point.SetAsTorsionFromString(s:String);
begin
Self.SetAsTorsionFromHash(SHA256StringHash(s));
end;

{ ******************************************************************************* }
procedure Fp4Point.SetAsRandomTorsionPoint;
var tmp: Fp4Point;
begin
tmp.SetCurveParams(CurveParams);
Infinity := false;
Randomize;
repeat
X.a.SetToRandom;
X.b.SetToZero;
SetPointValue(X);
if not Infinity then begin
                     if CurveParams.Family=cfKSS16 then _Mul_Fp4_FpPoint(CurveParams.htw, Self, tmp, csaffine)
                     // We use instead a faster approach proposed by Luis Dominiguez
                     else  _Fast_Mul_Htw_Fp4Point(Self,tmp);
                     Self := tmp;
                     Self.Z.SetToZero;
                     Self.z.a.a:=1;
                     end;
until not Infinity;
{_Mul_Fp4_FpPoint(CurveParams.Rtw, Self, tmp, csaffine);   // For Testing Torsion point Validity
if not tmp.Infinity then showmessage('problem! G2');}
end;

{ ******************************************************************************* }
procedure Fp4Point.SetPairingPointCoordinates(PtX, PtY: LInt);
begin
CurrentXp := PtX;
CurrentYp := PtY;
end;

{ ******************************************************************************* }
procedure Fp4Point.SetPointValue(Value: Fp4Int);
begin
_Evaluate_Fp4_Point(Value, Self);
end;

{ ******************************************************************************* }
procedure Fp4Point.SetToDefaultGenerator;
begin
Self.X:=CurveParams.BLS24TwistGeneratorX;
Self.Y:=CurveParams.BLS24TwistGeneratorY;
Self.Z.a.SetFormString('1');
Self.Z.b.SetFormString('0');
Self.SetCurveParams(CurveParams);
Self.Infinity:=false;
end;

{ ******************************************************************************* }
function Fp4Point.ToHexString: String;
begin
if Infinity then Result := 'Infinity'
else Result := '(' + X.ToHexString + ',' + Y.ToHexString + ')';
end;

{ ******************************************************************************* }
function Fp4Point.ToDecimalString: String;
begin
if Infinity then Result := 'Infinity'
else Result := '(' + X.ToDecimalString + ',' + Y.ToDecimalString + ')';
end;

{ ******************************************************************************* }
class operator Fp4Point.Add(Left, Right: Fp4Point): Fp4Point;
begin
_Add_Affine_Fp4_Point(Left, Right, Result);
end;

{ ******************************************************************************* }
class operator Fp4Point.Multiply(Left: LInt; Right { Q } : Fp4Point): Fp4Point;
begin
_Mul_Fp4_FpPoint(Left, Right, Result);
end;

{ ******************************************************************************* }
class operator Fp4Point.Subtract(Left, Right: Fp4Point): Fp4Point;
begin
_Sub_Fp4_Point(Left, Right, Result);
end;

{ ******************************************************************************* }
class operator Fp4Point.Negative(const Value: Fp4Point): Fp4Point;
begin
_Neg_Fp4_Point(Value, Result);
end;

{ ******************************************************************************* }
class operator Fp4Point.Equal(const Left, Right: Fp4Point): boolean;
begin
Result := _Are_Equals_Fp4Points(Left, Right);
end;

{ ******************************************************************************* }
class operator Fp4Point.NotEqual(const Left, Right: Fp4Point): boolean;
begin
Result := not _Are_Equals_Fp4Points(Left, Right);
end;

end.

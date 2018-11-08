unit ECCFp2;

interface
uses Fp2Arithmetic,System.SysUtils,LargeIntegers,Fp12Arithmetic,Fp6Arithmetic,HashFunctions, GeneralTypes,VCL.Dialogs;
 type

  Fp2Point=record     //  Definition of a Point on EC Over Fp2 (E(Fp2))
      public
      CurveParams: PtrCurveParams;
      X,Y,Z:Fp2Int;
      ComputeLigneValue:boolean;
      CurrentXp,CurrentYp:LInt;
      LineAtP:Fp12Int;
      Infinity:boolean;
      procedure SetCurveParams(Value: PtrCurveParams;ComputLmd:boolean=false);
      procedure SetPointValue(Value:Fp2Int);
      function ToDecimalString:String;
      function ToHexString:String;
      procedure InitFromStrings(Xa,Xb,Ya,Yb:String;Za:String='1';Zb:String='0');
      procedure SetAsRandomPoint;
      procedure SetAsRandomTorsionPoint;
      procedure SetToDefaultGenerator;
      procedure SetAsTorsionFromHash(hash:TBytes);
      procedure SetAsTorsionFromString(s:String);
      procedure SetPairingPointCoordinates(PtX,PtY:LInt);
      function CompressToArray:TBytes;
      function IsOnTheCurve:boolean;
      Procedure DeCompressFromArray(a:TBytes);
      function FrobeniusMap:Fp2Point;
      function FrobeniusMap2:Fp2Point;
      function FrobeniusMap3:Fp2Point;
      class operator Add(Left, Right: Fp2Point): Fp2Point;
      class operator Subtract(Left, Right: Fp2Point): Fp2Point;
      class operator Multiply(Left: LInt; Right: Fp2Point): Fp2Point;
      class operator Equal(const Left,Right: Fp2Point): Boolean;
      class operator NotEqual(const Left, Right: Fp2Point): Boolean;
      class operator Negative(const Value: Fp2Point): Fp2Point;
      end;

    procedure _Add_Jacobian_Fp2_Point(Left,Right:Fp2Point;var Result:Fp2Point;ComputeLignValue:boolean=false);
    procedure _Double_Jacobian_Fp2_Point(const Value:Fp2Point;var Result:Fp2Point;ComputeLignValue:boolean=false);
    procedure _Add_Projective_Fp2_Point(Left{Q},Right{T}:Fp2Point;var Result:Fp2Point;ComputeLignValue:boolean=false);
    procedure _Double_Projective_Fp2_Point(const Value:Fp2Point;var Result:Fp2Point;ComputeLigneValue:boolean=false);
    procedure _Add_Affine_Fp2_Point(Left,Right:Fp2Point;var Result:Fp2Point;ComputeLignValue:boolean=false);
    procedure _Double_Affine_Fp2_Point(const Value:Fp2Point;var Result:Fp2Point;ComputeLignValue:boolean=false);
    procedure _Double_Fp2_Point(Value:Fp2Point;var Result:Fp2Point;CoordSys:CoordinatesSystem=csAffine);
    procedure _Add_Fp2_Point(Left,Right:Fp2Point;var Result:Fp2Point;CoordSys:CoordinatesSystem=csAffine);
    procedure _Sub_Fp2_Point(Left,Right:Fp2Point;var Result:Fp2Point;CoordSys:CoordinatesSystem=csAffine);
    procedure _Evaluate_Fp2_Point(value:Fp2Int; var Result:Fp2Point);
    procedure _Neg_Fp2_Point(Value:Fp2Point;var Result:Fp2Point);
    function _Is_OnCurve_Fp2Point(value:Fp2Point):boolean;
    procedure _Frobenius_Map2(const Value:Fp2Point;var Result:Fp2Point);
    procedure _Frobenius_Map(const Value:Fp2Point;var Result:Fp2Point);
    procedure _Frobenius_Map3(const Value:Fp2Point;var Result:Fp2Point);
    procedure _Mul_Fp_Fp2Point(const Left:LInt;const Right:Fp2Point;var Result:Fp2Point;CoordSys:CoordinatesSystem=csAffine);
    procedure _Jacobian_To_Affine_Fp2Point(Value:Fp2Point;var Result:Fp2Point);
    procedure _Projective_To_Affine_Fp2Point(Value:Fp2Point;var Result:Fp2Point);
    function _Are_Equals_Fp2Points(left,right:Fp2Point):boolean;
    procedure _Fast_Mul_Htw_Fp2Point(const Left:LInt;const Qx0:Fp2Point;var Result:Fp2Point);    //Right Should be different than Right

 implementation

{*******************************************************************************}
            ///      Procedures for Elliptic Curves Arithmetic on Fp2
{*******************************************************************************}

        {**********   Add two Fp2 Points  using Affine coordinates *****************}
procedure _Add_Affine_Fp2_Point(Left,Right:Fp2Point;var Result:Fp2Point;ComputeLignValue:boolean=false);
var t: array[0..4] of Fp2Int;
begin
if Left.Infinity then Result:=Right
else if Right.Infinity then Result:=Left
else begin
     if _Compare_FP2(Left.X,Right.X)=0 then begin
                            if _Compare_FP2(Left.Y,Right.Y)=0 then _Double_Affine_Fp2_Point(Left,Result)
                            else begin
                                 Result.Z.a:=1;
                                 Result.Z.b:=0;
                                 Result.SetCurveParams(Right.CurveParams,ComputeLignValue);
                                 Result.Infinity:=true;
                                 end;
                            end
     else begin
          Result.SetCurveParams(Left.CurveParams,ComputeLignValue);
          _Sub_FP2(Right.Y,Left.Y,t[0]);
          _Sub_FP2(Right.X,Left.X,t[1]);
          _Inv_FP2(t[1],t[2]);
          _Mul_FP2(t[0],t[2],t[4]); // Lambda
          _Sqr_FP2(t[4],t[3]);
          if Result.ComputeLigneValue then begin
          if Left.CurveParams.TwistMode=twDType then begin
                                                           Result.LineAtP.a.a.a:=-Right.CurrentYp;
                                                           Result.LineAtP.a.a.b:=0;
                                                           Result.LineAtP.a.b.a.Data.i32[-1]:=0;
                                                           Result.LineAtP.a.b.b.Data.i32[-1]:=0;
                                                           Result.LineAtP.a.c.a.Data.i32[-1]:=0;
                                                           Result.LineAtP.a.c.b.Data.i32[-1]:=0;
                                                           _Mul_FP_FP2(Right.CurrentXp,t[4],Result.LineAtP.b.a);
                                                           _Mul_FP2(Left.x,t[4],t[0]);
                                                           _Sub_FP2(Left.y,t[0],Result.LineAtP.b.b);
                                                           Result.LineAtP.b.c.a.Data.i32[-1]:=0;
                                                           Result.LineAtP.b.c.b.Data.i32[-1]:=0;
                                                           end
          else begin
               Result.LineAtP.b.b.a:=-Right.CurrentYp;
               Result.LineAtP.b.b.b:=0;
               Result.LineAtP.b.a.a.Data.i32[-1]:=0;
               Result.LineAtP.b.a.b.Data.i32[-1]:=0;
               Result.LineAtP.a.c.a.Data.i32[-1]:=0;
               Result.LineAtP.a.c.b.Data.i32[-1]:=0;
               _Mul_FP_FP2(Right.CurrentXp,t[4],Result.LineAtP.a.b);
               _Mul_FP2(Left.x,t[4],t[0]);
               _Sub_FP2(Left.y,t[0],Result.LineAtP.a.a);
               Result.LineAtP.b.c.a.Data.i32[-1]:=0;
               Result.LineAtP.b.c.b.Data.i32[-1]:=0;
               end;
          end;
          _Sub_FP2(t[3],Left.X,Result.X);
          _Sub_FP2(Result.X,Right.X,Result.X);
          _Sub_FP2(Left.X,Result.X,t[2]);
          _Mul_FP2(t[4],t[2],Result.Y);
          _Sub_FP2(Result.Y,Left.Y,Result.Y);
          end;
     end;
end;

        {**********   Double an Fp2 Point using Affine coordinates *****************}
procedure _Double_Affine_Fp2_Point(const Value:Fp2Point;var Result:Fp2Point;ComputeLignValue:boolean=false);
var t: array[0..3] of Fp2Int;
    ValX,ValY:Fp2Int;
begin
if(Value.Infinity)or(Value.Y.a.Data.i16[-2]=0)and (Value.Y.b.Data.i16[-2]=0) then Result.Infinity:=true
else begin
     Result.Infinity:=false;
     Result.SetCurveParams(Value.CurveParams,ComputeLignValue);
     ValX:=Value.X;
     ValY:=Value.Y;
     _Sqr_FP2(Value.X,t[0]);
     _Mul_FP_FP2(LInt(3),t[0],t[1]);
     //_Add_FP2(t[1],Value.Curve.Atw,t[1]);
     _Mul_FP_FP2(LInt(2),Value.Y,t[3]);
     _Inv_FP2(t[3],t[3]);
     _Mul_FP2(t[1],t[3],t[2]);  // Lambda
     _Sqr_FP2(t[2],t[3]);
     _Sub_FP2(t[3],Value.X,t[3]);
     _Sub_FP2(t[3],Value.X,Result.X);
     _Sub_FP2(ValX,Result.x,t[3]);
     _Mul_FP2(t[2],t[3],Result.Y);
     _Sub_FP2(Result.y,ValY,Result.Y);
     if Result.ComputeLigneValue then begin
      if Value.CurveParams.TwistMode=twDType then begin
                                                        _Sub_LInt(Value.CurveParams.P,Value.CurrentYp,Result.LineAtP.a.a.a);
                                                        Result.LineAtP.a.a.b:=0;
                                                        Result.LineAtP.a.b.a.Data.i32[-1]:=0;
                                                        Result.LineAtP.a.b.b.Data.i32[-1]:=0;
                                                        Result.LineAtP.a.c.a.Data.i32[-1]:=0;
                                                        Result.LineAtP.a.c.b.Data.i32[-1]:=0;
                                                        _Mul_FP_FP2(Value.CurrentXp,t[2],Result.LineAtP.b.a);
                                                        _Mul_FP2(Valx,t[2],t[1]);
                                                        _Sub_FP2(Valy,t[1],t[1]);
                                                        Result.LineAtP.b.b:=t[1];
                                                        Result.LineAtP.b.c.a.Data.i32[-1]:=0;
                                                        Result.LineAtP.b.c.b.Data.i32[-1]:=0;
                                                        end
      else begin
           _Sub_LInt(Value.CurveParams.P,Value.CurrentYp,Result.LineAtP.b.b.a);
           Result.LineAtP.b.b.b:=0;
           Result.LineAtP.b.a.a.Data.i32[-1]:=0;
           Result.LineAtP.b.a.b.Data.i32[-1]:=0;
           Result.LineAtP.a.c.a.Data.i32[-1]:=0;
           Result.LineAtP.a.c.b.Data.i32[-1]:=0;
           _Mul_FP_FP2(Value.CurrentXp,t[2],Result.LineAtP.a.b);
           _Mul_FP2(Valx,t[2],t[1]);
           _Sub_FP2(Valy,t[1],t[1]);
           Result.LineAtP.a.a:=t[1];
           Result.LineAtP.b.c.a.Data.i32[-1]:=0;
           Result.LineAtP.b.c.b.Data.i32[-1]:=0;
           end;
        end;
     end
end;
        {**********   Add two Fp2 Points using Jacobian coordinates *****************}
        // (X/Z^2,Y/Z^3)        https://eprint.iacr.org/2010/354.pdf
procedure _Add_Jacobian_Fp2_Point(Left,Right:Fp2Point;var Result:Fp2Point;ComputeLignValue:boolean=false);
 var t: array[0..10] of Fp2Int;
     Zr2,Ly2:Fp2Int;
begin
if Left.Infinity then Result:=Right
else if Right.Infinity then Result:=Left
else begin
     if _Compare_FP2(Left.X,Right.X)=0 then begin
                            if _Compare_FP2(Left.Y,Right.Y)=0 then _Double_Jacobian_Fp2_Point(Left,Result)
                            else begin
                                 Result.Z.a:=1;
                                 Result.Z.b:=0;
                                 Result.SetCurveParams(Right.CurveParams,ComputeLignValue);
                                 Result.Infinity:=true;
                                 end;
                            end
     else begin
          Result.SetCurveParams(Left.CurveParams,ComputeLignValue);
          _Sqr_FP2(Right.Z,Zr2);
          _Mul_FP2(Left.X,Zr2,t[0]);
          _Add_FP2(Left.Y,Right.Z,t[2]);
          _Sqr_FP2(t[2],t[1]);
          _Sqr_FP2(Left.Y,Ly2);
          _Sub_FP2(t[1],Ly2,t[1]);
          _Sub_FP2(t[1],Zr2,t[2]);
          _Mul_FP2(t[2],Zr2,t[1]);
          _Sub_FP2(t[0],Right.X,t[2]);
          _Sqr_FP2(t[2],t[3]);
          _Mul_FP_FP2(LInt(4),t[3],t[4]);
          _Mul_FP2(t[4],t[2],t[5]);
          _Mul_FP_FP2(LInt(2),Right.Y,t[6]);
          _Sub_FP2(t[1],t[6],t[6]);
          _Mul_FP2(t[6],Left.X,t[9]);
          _Mul_FP2(Right.X,t[4],t[7]);
          _Mul_FP_FP2(LInt(2),t[7],t[10]);
          _Sqr_FP2(t[6],Result.X);
          _Sub_FP2(Result.X,t[5],Result.x);
          _Sub_FP2(Result.X,t[10],Result.X);
          _Add_FP2(Right.Z,t[2],Result.Z);
          _Sqr_FP2(Result.Z,t[10]);
          _Sub_FP2(t[10],Zr2,Result.Z);
          _Sub_FP2(Result.Z,t[3],Result.Z);
          _Add_FP2(Left.Y,Result.Z,t[10]);
          _Sub_FP2(t[7],Result.X,t[7]);
          _Mul_FP2(t[7],t[6],t[8]);
          _Mul_FP2(Right.Y,t[5],t[0]);
          _Shl_LInt(t[0].a,1);
          _Shl_LInt(t[0].b,1);
          _Mod_LInt(t[0].a,Left.CurveParams.P,t[0].a);
          _Mod_LInt(t[0].b,Left.CurveParams.P,t[0].b);
          _Sub_FP2(t[8],t[0],Result.Y);
          if Result.ComputeLigneValue then begin
          if Left.CurveParams.TwistMode=twDType then begin
                                                           _Sqr_FP2(t[10],t[3]);
                                                           _Sub_FP2(t[3],Ly2,t[3]);
                                                           _Sqr_FP2(Result.Z,t[4]);
                                                           _Sub_FP2(t[3],t[4],t[10]);
                                                           _Shl_LInt(t[9].a,1);
                                                           _Shl_LInt(t[9].b,1);
                                                           _Mod_LInt(t[9].a,Right.CurveParams.P,t[9].a);
                                                           _Mod_LInt(t[9].b,Right.CurveParams.P,t[9].b);
                                                           _Sub_FP2(t[9],t[10],t[9]);
                                                           _Mul_FP_FP2(Right.CurrentYp,Result.Z,t[10]);
                                                           _Shl_LInt(t[10].a,1);
                                                           _Shl_LInt(t[10].b,1);
                                                           _Mod_LInt(t[10].a,Right.CurveParams.P,t[10].a);
                                                           _Mod_LInt(t[10].b,Right.CurveParams.P,t[10].b);
                                                           _Sub_LInt(Right.CurveParams.P,t[6].a,t[6].a);
                                                           _Sub_LInt(Right.CurveParams.P,t[6].b,t[6].b);
                                                           _Mul_FP_FP2(Right.CurrentXp,t[6],t[1]);
                                                           _Shl_LInt(t[1].a,1);
                                                           _Shl_LInt(t[1].b,1);
                                                           _Mod_LInt(t[1].a,Right.CurveParams.P,t[1].a);
                                                           _Mod_LInt(t[1].b,Right.CurveParams.P,t[1].b);
                                                           Result.LineAtP.a.a:=t[10];
                                                           Result.LineAtP.a.b.a.Data.i32[-1]:=0;
                                                           Result.LineAtP.a.b.b.Data.i32[-1]:=0;
                                                           Result.LineAtP.a.c.a.Data.i32[-1]:=0;
                                                           Result.LineAtP.a.c.b.Data.i32[-1]:=0;
                                                           Result.LineAtP.b.a:=t[1];
                                                           Result.LineAtP.b.b:=t[9];
                                                           Result.LineAtP.b.c.a.Data.i32[-1]:=0;
                                                           Result.LineAtP.b.c.b.Data.i32[-1]:=0;
                                                           end
          else begin
               _Sqr_FP2(t[10],t[3]);
               _Sub_FP2(t[3],Ly2,t[3]);
               _Sqr_FP2(Result.Z,t[4]);
               _Sub_FP2(t[3],t[4],t[10]);
               _Shl_LInt(t[9].a,1);
               _Shl_LInt(t[9].b,1);
               _Mod_LInt(t[9].a,Right.CurveParams.P,t[9].a);
               _Mod_LInt(t[9].b,Right.CurveParams.P,t[9].b);
               _Sub_FP2(t[9],t[10],t[9]);
               _Mul_FP_FP2(Right.CurrentYp,Result.Z,t[10]);
               _Shl_LInt(t[10].a,1);
               _Shl_LInt(t[10].b,1);
               _Mod_LInt(t[10].a,Right.CurveParams.P,t[10].a);
               _Mod_LInt(t[10].b,Right.CurveParams.P,t[10].b);
               _Sub_LInt(Right.CurveParams.P,t[6].a,t[6].a);
               _Sub_LInt(Right.CurveParams.P,t[6].b,t[6].b);
               _Mul_FP_FP2(Right.CurrentXp,t[6],t[1]);
               _Shl_LInt(t[1].a,1);
               _Shl_LInt(t[1].b,1);
               _Mod_LInt(t[1].a,Right.CurveParams.P,t[1].a);
               _Mod_LInt(t[1].b,Right.CurveParams.P,t[1].b);
               Result.LineAtP.b.b:=t[10];
               Result.LineAtP.b.a.a.Data.i32[-1]:=0;
               Result.LineAtP.b.a.b.Data.i32[-1]:=0;
               Result.LineAtP.a.c.a.Data.i32[-1]:=0;
               Result.LineAtP.a.c.b.Data.i32[-1]:=0;
               Result.LineAtP.a.b:=t[1];
               Result.LineAtP.a.a:=t[9];
               Result.LineAtP.b.c.a.Data.i32[-1]:=0;
               Result.LineAtP.b.c.b.Data.i32[-1]:=0;
               end;
            end;
          end;
     end;
end;
        {**********   Double an Fp2 Point using Jacobian coordinates *****************}
        // (X/Z^2,Y/Z^3)        https://eprint.iacr.org/2010/354.pdf
procedure _Double_Jacobian_Fp2_Point(const Value:Fp2Point;var Result:Fp2Point;ComputeLignValue:boolean=false);
var t: array[0..7] of Fp2Int;
    ZQ2:Fp2Int;
begin
if(Value.Infinity)or((Value.Y.a.Data.i16[-2]=0) and (Value.Y.b.Data.i16[-2]=0)) then Result.Infinity:=true
else begin
     Result.Infinity:=false;
     Result.SetCurveParams(Value.CurveParams,ComputeLignValue);
     _Sqr_FP2(Value.Z,ZQ2);
     _Sqr_FP2(Value.X,t[0]);
     _Sqr_FP2(Value.Y,t[1]);
     _Sqr_FP2(t[1],t[2]);
     _Add_FP2(t[1],Value.X,t[3]);
     _Sqr_FP2(t[3],t[4]);
     _Sub_FP2(t[4],t[0],t[3]);
     _Sub_FP2(t[3],t[2],t[3]);
     _Shl_LInt(t[3].a,1);
     _Shl_LInt(t[3].b,1);
     _Mod_LInt(t[3].a,Value.CurveParams.P,t[3].a);
     _Mod_LInt(t[3].b,Value.CurveParams.P,t[3].b);
     _Mul_FP_FP2(LInt(3),t[0],t[4]);
     _Add_FP2(Value.X,t[4],t[6]);
     _Sqr_FP2(t[4],t[5]);
     _Mul_FP_FP2(LInt(2),t[3],Result.X);
     _Sub_FP2(t[5],Result.X,Result.X);
     _Add_FP2(Value.Y,Value.Z,t[7]);
     _Sqr_FP2(t[7],Result.Z);
     _Sub_FP2(Result.Z,t[1],Result.Z);
     _Sub_FP2(Result.Z,ZQ2,Result.Z);
     _Sub_FP2(t[3],Result.X,Result.Y);
     _Mul_FP2(Result.Y,t[4],t[7]);
     _Mul_FP_FP2(LInt(8),t[2],Result.Y);
     _Sub_FP2(t[7],Result.Y,Result.Y);
     if Result.ComputeLigneValue then begin
     if Value.CurveParams.TwistMode=twDType then begin
                                                       _Mul_FP2(t[4],ZQ2,t[3]);
                                                       _Mul_FP_FP2(LInt(2),t[3],t[7]);
                                                       _Sub_LInt(Value.CurveParams.P,t[7].a,t[7].a);
                                                       _Sub_LInt(Value.CurveParams.P,t[7].b,t[7].b);
                                                       _Mul_FP_FP2(Value.CurrentXp,t[7],t[3]);
                                                       _Sqr_FP2(t[6],t[7]);
                                                       _Mul_FP_FP2(LInt(4),t[1],t[2]);
                                                       _Sub_FP2(t[7],t[0],t[7]);
                                                       _Sub_FP2(t[7],t[5],t[7]);
                                                       _Sub_FP2(t[7],t[2],t[6]);
                                                       _Mul_FP2(Result.Z,ZQ2,t[2]);
                                                       _Mul_FP_FP2(LInt(2),t[2],t[7]);
                                                       _Mul_FP_FP2(Value.CurrentYp,t[7],t[0]);
                                                       Result.LineAtP.a.a:=t[0];
                                                       Result.LineAtP.a.b.a.Data.i32[-1]:=0;
                                                       Result.LineAtP.a.b.b.Data.i32[-1]:=0;
                                                       Result.LineAtP.a.c.a.Data.i32[-1]:=0;
                                                       Result.LineAtP.a.c.b.Data.i32[-1]:=0;
                                                       Result.LineAtP.b.a:=t[3];
                                                       Result.LineAtP.b.b:=t[6];
                                                       Result.LineAtP.b.c.a.Data.i32[-1]:=0;
                                                       Result.LineAtP.b.c.b.Data.i32[-1]:=0;
                                                       end
     else begin
          _Mul_FP2(t[4],ZQ2,t[3]);
          _Mul_FP_FP2(LInt(2),t[3],t[7]);
          _Sub_LInt(Value.CurveParams.P,t[7].a,t[7].a);
          _Sub_LInt(Value.CurveParams.P,t[7].b,t[7].b);
          _Mul_FP_FP2(Value.CurrentXp,t[7],t[3]);
          _Sqr_FP2(t[6],t[7]);
          _Mul_FP_FP2(LInt(4),t[1],t[2]);
          _Sub_FP2(t[7],t[0],t[7]);
          _Sub_FP2(t[7],t[5],t[7]);
          _Sub_FP2(t[7],t[2],t[6]);
          _Mul_FP2(Result.Z,ZQ2,t[2]);
          _Mul_FP_FP2(LInt(2),t[2],t[7]);
          _Mul_FP_FP2(Value.CurrentYp,t[7],t[0]);
          Result.LineAtP.b.b:=t[0];
          Result.LineAtP.b.a.a.Data.i32[-1]:=0;
          Result.LineAtP.b.a.b.Data.i32[-1]:=0;
          Result.LineAtP.a.c.a.Data.i32[-1]:=0;
          Result.LineAtP.a.c.b.Data.i32[-1]:=0;
          Result.LineAtP.a.b:=t[3];
          Result.LineAtP.a.a:=t[6];
          Result.LineAtP.b.c.a.Data.i32[-1]:=0;
          Result.LineAtP.b.c.b.Data.i32[-1]:=0;
          end;
        end;
     end
end;

        {**********   Double an Fp2 Point (Projective coordinates (X/Z,Y/Z))*****************}
        // Costello et al's   https://cryptojedi.org/papers/edate-20100614.pdf
procedure _Double_Projective_Fp2_Point(const Value:Fp2Point;var Result:Fp2Point;ComputeLigneValue:boolean=false);
var tmp,A,B,C,D,E,F,G:Fp2Int;
begin
if(Value.Infinity)or((Value.Y.a.Data.i16[-2]=0)and (Value.Y.b.Data.i16[-2]=0))or(_Is_FP2_Null(Value.Z)) then Result.Infinity:=true
else begin
     Result.Infinity:=false;
     Result.SetCurveParams(Value.CurveParams,ComputeLigneValue);
     _Sqr_FP2(Value.X,A);
     _Sqr_FP2(Value.Y,B);
     _Sqr_FP2(Value.Z,C);
     _Mul_FP2(C,Value.CurveParams.Btw,tmp);
     _Nomod_Add_FP2(tmp,tmp,D);
     _Add_FP2(tmp,D,D);
     _Nomod_Add_FP2(Value.X,Value.Y,tmp);
     _Sqr_FP2(tmp,E);
     _Nomod_Sub_FP2(E,A,E);
     _Sub_FP2(E,B,E);
     _Nomod_Add_FP2(Value.Y,Value.Z,tmp);
     _Sqr_FP2(tmp,F);
     _Nomod_Sub_FP2(F,B,F);
     _Sub_FP2(F,C,F);
     _Nomod_Add_FP2(D,D,G);
     _Add_FP2(D,G,G);
     _Nomod_Sub_FP2(B,G,tmp);
     _Mul_FP2(tmp,E,Result.X);
     _Add_FP2(B,G,tmp);
     _Sqr_FP2(tmp,Result.Y);
     _Sqr_FP2(D,tmp);
     _Shl_LInt(tmp.a,2);
     _Shl_LInt(tmp.b,2);
     _Sub_FP2(Result.Y,tmp,Result.Y);
     _Nomod_Add_FP2(tmp,tmp,tmp);
     _Sub_FP2(Result.Y,tmp,Result.Y);
     _Mul_FP2(B,F,Result.Z);
     _Shl_LInt(Result.Z.a,2);
     _Shl_LInt(Result.Z.b,2);
     _Mod_LInt(Result.Z.a,Value.CurveParams.P,Result.Z.a);
     _Mod_LInt(Result.Z.b,Value.CurveParams.P,Result.Z.b);
     if Result.ComputeLigneValue then begin
     if Value.CurveParams.TwistMode=twDType then begin
                                                       _Mul_FP_FP2(Value.CurrentYp,F,tmp);
                                                       _Neg_FP2(tmp,Result.LineAtP.a.a);     //l11
                                                       Result.LineAtP.a.b.a.Data.i32[-1]:=0;
                                                       Result.LineAtP.a.b.b.Data.i32[-1]:=0;
                                                       Result.LineAtP.a.c.a.Data.i32[-1]:=0;
                                                       Result.LineAtP.a.c.b.Data.i32[-1]:=0;
                                                        _Nomod_Add_FP2(A,A,tmp);
                                                       _Add_FP2(A,tmp,tmp);
                                                       _Mul_FP_FP2(Value.CurrentXp,tmp,Result.LineAtP.b.a);    //l02
                                                        _Sub_FP2(D,B, Result.LineAtP.b.b);     //l00
                                                       Result.LineAtP.b.c.a.Data.i32[-1]:=0;
                                                       Result.LineAtP.b.c.b.Data.i32[-1]:=0;
                                                       end
     else begin
          _Mul_FP_FP2(Value.CurrentYp,F,tmp);
          _Neg_FP2(tmp,Result.LineAtP.b.b);     //l11
          Result.LineAtP.b.a.a.Data.i32[-1]:=0;
          Result.LineAtP.b.a.b.Data.i32[-1]:=0;
          Result.LineAtP.a.c.a.Data.i32[-1]:=0;
          Result.LineAtP.a.c.b.Data.i32[-1]:=0;
         _Nomod_Add_FP2(A,A,tmp);
         _Add_FP2(A,tmp,tmp);
         _Mul_FP_FP2(Value.CurrentXp,tmp,Result.LineAtP.a.b);    //l02
         _Sub_FP2(D,B, Result.LineAtP.a.a);     //l00
         Result.LineAtP.b.c.a.Data.i32[-1]:=0;
         Result.LineAtP.b.c.b.Data.i32[-1]:=0;
         end;
      end;
   end
end;

        {**********   Add two Fp2 Points (Projective coordinates)*****************}
        // Mixed addition (Left.Z is always equal to 1)
        // Costello et al's   https://cryptojedi.org/papers/edate-20100614.pdf
procedure _Add_Projective_Fp2_Point(Left{Q},Right{T}:Fp2Point;var Result:Fp2Point;ComputeLignValue:boolean=false);
var A,B,C,D,tmp1,tmp2:Fp2Int;
begin
if Left.Infinity then Result:=Right
else if Right.Infinity then Result:=Left
else begin
     if _Compare_FP2(Left.X,Right.X)=0 then begin
                            if _Compare_FP2(Left.Y,Right.Y)=0 then _Double_Projective_Fp2_Point(Left,Result)
                            else begin
                                 Result.Infinity:=true;
                                 Result.Z.a:=1;
                                 Result.Z.b:=0;
                                 Result.SetCurveParams(Right.CurveParams,ComputeLignValue);
                                 end;
                            end
     else begin
          Result.SetCurveParams(Left.CurveParams,ComputeLignValue);
          _Mul_FP2(Right.Z,Left.X,A);
          _Sub_FP2(Right.X,A,A);
          _Mul_FP2(Right.Z,Left.Y,B);
          _Sub_FP2(Right.Y,B,B);
          _Sqr_FP2(A,tmp1);
          _Mul_FP2(Right.X,tmp1,tmp2); // tmp2 used later as Right.X
          _Mul_FP2(tmp1,A,C);
          _Mul_FP2(B,B,tmp1);
          _Mul_FP2(tmp1,Right.Z,D);
          _Nomod_Add_FP2(D,C,D);
          _Add_FP2(tmp2,tmp2,tmp1);
          _Sub_FP2(D,tmp1,D);
          _Sub_FP2(tmp2,D,tmp2);
          _Mul_FP2(B,tmp2,Result.Y);
          _Mul_FP2(Right.Y,C,tmp2);
          _Sub_FP2(Result.Y,tmp2,Result.Y);
          _Mul_FP2(A,D,Result.X);
          _Mul_FP2(Right.Z,C,Result.Z);
          if Result.ComputeLigneValue then begin
          if  Left.CurveParams.TwistMode=twDType then begin
                                                           _Mul_FP_FP2(Right.CurrentYp,A,Result.LineAtP.a.a);      //L0,1
                                                           Result.LineAtP.a.b.a.Data.i32[-1]:=0;
                                                           Result.LineAtP.a.b.b.Data.i32[-1]:=0;
                                                           Result.LineAtP.a.c.a.Data.i32[-1]:=0;
                                                           Result.LineAtP.a.c.b.Data.i32[-1]:=0;
                                                           _Mul_FP_FP2(Right.CurrentXp,B,tmp1);
                                                           _Neg_FP2(tmp1,Result.LineAtP.b.a);
                                                           _Mul_FP2(B,Left.X,tmp1);
                                                           _Mul_FP2(A,Left.Y,tmp2);
                                                           _Sub_FP2(tmp1,tmp2,Result.LineAtP.b.b);
                                                           Result.LineAtP.b.c.a.Data.i32[-1]:=0;
                                                           Result.LineAtP.b.c.b.Data.i32[-1]:=0;
                                                           end
          else begin
               _Mul_FP_FP2(Right.CurrentYp,A,Result.LineAtP.b.b);      //L0,1
               Result.LineAtP.b.a.a.Data.i32[-1]:=0;
               Result.LineAtP.b.a.b.Data.i32[-1]:=0;
               Result.LineAtP.a.c.a.Data.i32[-1]:=0;
               Result.LineAtP.a.c.b.Data.i32[-1]:=0;
               _Mul_FP_FP2(Right.CurrentXp,B,tmp1);
               _Neg_FP2(tmp1,Result.LineAtP.a.b);
               _Mul_FP2(B,Left.X,tmp1);
               _Mul_FP2(A,Left.Y,tmp2);
               _Sub_FP2(tmp1,tmp2,Result.LineAtP.a.a);
               Result.LineAtP.b.c.a.Data.i32[-1]:=0;
               Result.LineAtP.b.c.b.Data.i32[-1]:=0;
               end;
            end;
          end;
     end;
end;

        {**********   Add two Fp2 Points *****************}
procedure _Add_Fp2_Point(Left,Right:Fp2Point;var Result:Fp2Point;CoordSys:CoordinatesSystem=csAffine);
begin
case CoordSys of
  csAffine:_Add_Affine_Fp2_Point(Left,Right,Result);
  csProjective:_Add_Projective_Fp2_Point(Left,Right,Result);
  csJacobian:_Add_Jacobian_Fp2_Point(Left,Right,Result);
end;
end;

        {**********   Double two Fp2 Points *****************}
procedure _Double_Fp2_Point(Value:Fp2Point;var Result:Fp2Point;CoordSys:CoordinatesSystem=csAffine);
begin
case CoordSys of
  csAffine:_Double_Affine_Fp2_Point(Value,Result);
  csProjective:_Double_Projective_Fp2_Point(Value,Result);
  csJacobian:_Double_Jacobian_Fp2_Point(Value,Result);
end;
end;

        {**********   Substract two Fp2 Points *****************}
procedure _Sub_Fp2_Point(Left,Right:Fp2Point;var Result:Fp2Point;CoordSys:CoordinatesSystem=csAffine);
var tmp:Fp2Point;
begin
_neg_fp2_point(Right,tmp);
case CoordSys of
  csAffine:_Add_Affine_Fp2_Point(Left,tmp,Result);
  csProjective:_Add_Projective_Fp2_Point(Left,tmp,Result);
  csJacobian:_Add_Jacobian_Fp2_Point(Left,tmp,Result);
end;
end;

       {**********   Convert Jacobian/Affine for an Fp2 Point *****}
        //(X/Z^2,Y/Z^3)
procedure _Jacobian_To_Affine_Fp2Point(Value:Fp2Point;var Result:Fp2Point);
var t1,t2:Fp2Int;
begin
if (_Is_FP2_Null(Value.Z)) then Result.Infinity:=true
else begin
     t1:=Value.z.Inverse;
     t2:=t1.Sqr;
     Result.X:=(Value.X*t2);
     Result.y:=(Value.y*(t1*t2));
     Result.z.a:=1;
     Result.z.b:=0;
     end;
end;

       {**********   Convert Projective/Affine for an Fp2 Point *****}
       // (X/Z,Y/Z)
procedure _Projective_To_Affine_Fp2Point(Value:Fp2Point;var Result:Fp2Point);
var t1,t2:Fp2Int;
begin
if (_Is_FP2_Null(Value.Z)) then Result.Infinity:=true
else begin
     t1:=Value.z.Inverse;
     t2:=t1.Sqr;
     Result.X:=(Value.X*t1);
     Result.y:=(Value.y*t1);
     Result.z.a:=1;
     Result.z.b:=0;
     end;
end;

       {**********   Normalize coordinates to Affine for an Fp2 Point *****}
procedure _Normalize_PF2_Point(value:FP2Point;var Result:FP2Point;CoordSys:CoordinatesSystem=csAffine);
begin
case  CoordSys of
  csAffine:result:=value;
  csJacobian:_Jacobian_To_Affine_Fp2Point(value,Result);
  csProjective:_Projective_To_Affine_Fp2Point(Value,Result);
end;
end;

      { **********   Multiply the Cofactor with an Fp2 Point ******* }
      // Fast Multiplication of A point by the Twist Cofactor (Generated By  Luis Dominguez)
      // For BLS12 Curves (Not for BN curves)
      // https://eprint.iacr.org/2008/530.pdf
procedure _Fast_Mul_Htw_Fp2Point(const Left:LInt;const Qx0:Fp2Point;var Result:Fp2Point);    //Right Should be different than Result
var xA,xB,xC,t0,Qx0_,Qx1_,Qx1,Qx2_,Qx2,Qx3_,Qx3:Fp2Point;
begin
Result.CurveParams :=Qx0.CurveParams;
if Qx0.Infinity then begin
                     Result.Infinity:=true;
                     exit;
                     end;
_Neg_Fp2_Point(Qx0,Qx0_);
_Mul_Fp_Fp2Point(Left,Qx0,Qx1);
_Neg_Fp2_Point(Qx1,Qx1_);
_Mul_Fp_Fp2Point(Left,Qx1,Qx2);
_Neg_Fp2_Point(Qx2,Qx2_);
_Mul_Fp_Fp2Point(Left,Qx2,Qx3);
_Neg_Fp2_Point(Qx3,Qx3_);
xA:=Qx0;
_Add_Affine_Fp2_Point(xA,Qx0,t0);
_Frobenius_Map2(Qx1,xB);
_Add_Affine_Fp2_Point(t0,xB,t0);
_Add_Affine_Fp2_Point(t0,t0,t0);
xB:=Qx1_;
xC:=Qx2_;
_Add_Affine_Fp2_Point(xB,xC,xB);
xC:=Qx3;
_Add_Affine_Fp2_Point(xB,xC,xB);
_Frobenius_Map(Qx0,xC);
_Add_Affine_Fp2_Point(xB,xC,xB);
_Frobenius_Map(Qx1_,xC);
_Add_Affine_Fp2_Point(xB,xC,xB);
_Frobenius_Map(Qx2_,xC);
_Add_Affine_Fp2_Point(xB,xC,xB);
_Frobenius_Map(Qx3,xC);
_Add_Affine_Fp2_Point(xB,xC,xB);
_Frobenius_Map2(Qx0_,xC);
_Add_Affine_Fp2_Point(xB,xC,xB);
_Frobenius_Map2(Qx2_,xC);
_Add_Affine_Fp2_Point(xB,xC,xB);
_Add_Affine_Fp2_Point(t0,xB,Result);
Result.Z.a:=1;
Result.Z.b:=0;
Result.Z.Field:=Qx0.CurveParams.FieldParam;
end;

       {**********   Multiply a scalar with an Fp2 Point *******}
procedure _Mul_Fp_Fp2Point(const Left:LInt;const Right:Fp2Point;var Result:Fp2Point; CoordSys:CoordinatesSystem=csAffine);    //Right Should be different than Right
var i:integer;
begin
Result.CurveParams:=Right.CurveParams;
if (Right.Infinity)or (Left.Data.i16[-2]=0) then Result.Infinity:=true
else begin
     Result.Infinity:=true;
     if Left=2 then begin
                    case CoordSys of
                    csAffine:_Double_Affine_Fp2_Point(Right,Result);
                    csProjective:_Double_Projective_Fp2_Point(Right,Result);
                    csJacobian:_Double_Jacobian_Fp2_Point(Right,Result);
                    end;
                    end
     else
     for i:=Left.BitLength-1 downto 0 do begin
                                            case CoordSys of
                                              csAffine:_Double_Affine_Fp2_Point(Result,Result);
                                              csProjective:_Double_Projective_Fp2_Point(Result,Result);
                                              csJacobian:_Double_Jacobian_Fp2_Point(Result,Result);
                                            end;
                                         if _Is_BitSet_At(Left,i) then begin
                                                                  case CoordSys of
                                                                    csAffine:_Add_Affine_Fp2_Point(Right,Result,Result);
                                                                    csProjective:_Add_Projective_Fp2_Point(Right,Result,Result);
                                                                    csJacobian:_Add_Jacobian_Fp2_Point(Right,Result,Result);
                                                                  end;
                                                                 end;
                                         end;
     _Normalize_PF2_Point(Result,Result,CoordSys);
     if _IsNeg(left) then _Neg_Fp2_Point(Result,Result);
     end;
end;

       {**********   Multiply a scalar with an Fp Point *******}
procedure _Mul_NAF_Fp2_FpPoint(const Left:LInt;const Right:Fp2Point;var Result:Fp2Point;CoordSys:CoordinatesSystem=csAffine);
var i:integer;
    loop:LIntArrayForm;
    NegLeft:Fp2Point;
begin
Result.CurveParams:=Right.CurveParams;
_Neg_Fp2_Point(Right,NegLeft);
if (Right.Infinity)or (Left.Data.i16[-2]=0) then Result.Infinity:=true
else begin
     Result.Infinity:=true;
     if Left=2 then _Double_Jacobian_Fp2_Point(Right,Result)
     else begin
          Loop:=Left.ToNafArray;
          for i:=Length(Loop)-1 downto 0 do begin
                                            _Double_Fp2_Point(Result,Result,CoordSys);
                                            if Loop[i]=1 then
                                                _Add_Fp2_Point(Right,Result,Result,CoordSys);
                                            if Loop[i]=-1 then
                                                _Add_Fp2_Point(NegLeft,Result,Result,CoordSys);
                                            end;
          _Normalize_PF2_Point(Result,Result,CoordSys);
          end;
     end;
end;

       {**********   Compute Frobenius Map for an Fp2 Point *******}
procedure _Frobenius_Map(const Value:Fp2Point;var Result:Fp2Point);
begin
if Value.Infinity then Result.Infinity:=true
else begin
     Result.SetCurveParams (Value.CurveParams);
     _Mul_Fp2(Value.CurveParams.FrobeniusMapConstX[0],Value.X.Conjugate,Result.X);
     _Mul_Fp2(Value.CurveParams.FrobeniusMapConstY[0],Value.Y.Conjugate,Result.Y);
     _HCopy_LInt(Value.Z.a,Result.Z.a);
     _HCopy_LInt(Value.Z.b,Result.Z.b);
     end;
end;

       {**********   Compute Frobenius2 Map for an Fp2 Point *******}
procedure _Frobenius_Map2(const Value:Fp2Point;var Result:Fp2Point);
begin
if Value.Infinity then Result.Infinity:=true
else begin
     Result.SetCurveParams (Value.CurveParams);
     _Mul_Fp2(Value.CurveParams.FrobeniusMapConstX[1],Value.X,Result.X);
     _Mul_Fp2(Value.CurveParams.FrobeniusMapConstY[1],Value.Y,Result.Y);
     _HCopy_LInt(Value.Z.a,Result.Z.a);
     _HCopy_LInt(Value.Z.b,Result.Z.b);
     end;
end;

       {**********   Compute Frobenius3 Map for an Fp2 Point *******}
procedure _Frobenius_Map3(const Value:Fp2Point;var Result:Fp2Point);
begin
if Value.Infinity then Result.Infinity:=true
else begin
     Result.SetCurveParams(Value.CurveParams);
     _Mul_Fp2(Value.CurveParams.FrobeniusMapConstX[2],Value.X.Conjugate,Result.X);
     _Mul_Fp2(Value.CurveParams.FrobeniusMapConstY[2],Value.Y.Conjugate,Result.Y);
     _HCopy_LInt(Value.Z.a,Result.Z.a);
     _HCopy_LInt(Value.Z.b,Result.Z.b);
     end;
end;
   {**********   Test if an Fp Point is on the Curve ***********}
function _Is_OnCurve_Fp2Point(value:Fp2Point):boolean;
var tmp,tmp1:Fp2Int;
begin
if value.Infinity then Result:=true
else begin
     _Sqr_FP2(Value.X,tmp);
     _Add_FP2(tmp,Value.CurveParams.Atw,tmp);
     _Mul_FP2(tmp,Value.X,tmp1);
     _Add_FP2(tmp1,Value.CurveParams.Btw,tmp);
     _Sqr_FP2(Value.Y,tmp1);
    Result:=(tmp.IsASquare) and (tmp1=tmp);
    end;
end;

        {**********   Negation of an Fp2 Point *************}
procedure _Neg_Fp2_Point(Value:Fp2Point;var Result:Fp2Point);
begin
  Result:=Value;
  _Neg_FP2(Result.Y,Result.Y);
end;

       {**********   Test if two Fp2 Points are equal *******}
function _Are_Equals_Fp2Points(left,right:Fp2Point):boolean;
begin
Result:=(_Compare_FP2(Left.X,Right.X)=0)and(_Compare_FP2(Left.Y,Right.Y)=0);
end;

      {**********   Evaluate an Fp2 Point from X value ******}
procedure _Evaluate_Fp2_Point(value:Fp2Int; var Result:Fp2Point);
var tmp,tmp1:Fp2Int;
begin
Result.X:=Value;
Result.Infinity:=false;
_Sqr_FP2(Result.X,tmp);
_Add_FP2(tmp,Result.CurveParams.Atw,tmp);
_Mul_FP2(tmp,Result.X,tmp1);
_Add_FP2(tmp1,Result.CurveParams.Btw,tmp);
if tmp.IsASquare then begin
                      _Sqrt_FP2(tmp,Result.Y);
                      _Neg_FP2(Result.Y,tmp1);
                      if _Compare_FP2(Result.Y,tmp1)=1 then Result.Y:=tmp1;
                      Result.Z.a:=1;
                      Result.Z.b:=0;
                      end
else Result.Infinity:=true;
end;

{*******************************************************************************}
      ///  Definitions of an Fp2 Point operators and functions
{*******************************************************************************}

{procedure Fp2Point.SetCurve(value: TCurve);
begin
_Curve:=Value;
if (Value<>nil)and(Value.CurveParams<>nil) then SetCurveParams(Value.CurveParams);
end;}
{*******************************************************************************}
procedure Fp2Point.SetCurveParams(Value: PtrCurveParams;ComputLmd:boolean=false);
begin
CurveParams:=Value;
if CurveParams.TowerParam<>nil then LineAtP.SetTowerParams (Value.TowerParam);
X.SetFieldParams(Value.FieldParam);
Y.SetFieldParams(Value.FieldParam);
Z.SetFieldParams(Value.FieldParam);
ComputeLigneValue:=ComputLmd;
Infinity:=false;
end;
{*******************************************************************************}
function Fp2Point.CompressToArray: TBytes;
var L:integer;
begin
L:=X.Field.p.Data.i32[-1]*4*2;
Setlength(Result,L);
Move(X.a.Data.i8[0],Result[0],Length(Result)shr 1);
Move(X.b.Data.i8[0],Result[L shr 1],Length(Result) shr 1);
end;
{*******************************************************************************}
procedure Fp2Point.DeCompressFromArray(a: TBytes);
begin
Move(a[0],X.a.Data.i8[0],Length(a) shr 1);
Move(a[Length(a) shr 1],X.b.Data.i8[0],Length(a) shr 1);
X.a.Data.i16[-2]:=Length(a) shr 3;
X.a.Data.i16[-1]:=0;
X.b.Data.i16[-2]:=Length(a) shr 3;
X.b.Data.i16[-1]:=0;
SetPointValue(X);
end;
{*******************************************************************************}
function Fp2Point.FrobeniusMap: Fp2Point;
begin
_Frobenius_Map(Self,Result);
end;
{*******************************************************************************}
function Fp2Point.FrobeniusMap2: Fp2Point;
begin
_Frobenius_Map2(Self,Result);
end;
{*******************************************************************************}
function Fp2Point.FrobeniusMap3: Fp2Point;
begin
_Frobenius_Map3(Self,Result);
end;
{*******************************************************************************}
procedure Fp2Point.InitFromStrings(Xa,Xb,Ya,Yb:String;Za:String='1';Zb:String='0');
begin
X.a:=Xa;
X.b:=Xb;
Y.a:=Ya;
Y.b:=Yb;
Z.a:=Za;
Z.b:=Zb;
Infinity:=false;
end;
{*******************************************************************************}
function Fp2Point.IsOnTheCurve: boolean;
begin
Result:=_Is_OnCurve_Fp2Point(Self)
end;
{*******************************************************************************}
procedure Fp2Point.SetAsRandomPoint;
begin
Infinity:=false;
Randomize;
repeat
  GetRandomLIntLowerThan(X.a,CurveParams.P);
  GetRandomLIntLowerThan(X.b,CurveParams.P);
  SetPointValue(X);
  until not Infinity;
end;
{*******************************************************************************}
procedure Fp2Point.SetAsTorsionFromHash(hash: TBytes);
var i,j:integer;
    h1,h2:TBytes;
    tmp:Fp2Point;
begin
j:=0;
repeat
  h1:=SHA256BytesHash(hash);
  h1[0]:=j;
  h2:=SHA256BytesHash(h1);
  Infinity:=false;
  X.a.Data.i32[length(h1) div 4]:=0;
  for i:=0 to length(h1)-1 do X.a.Data.i8[i]:=h1[i];
  X.a.Data.i16[-1]:=0;
  X.a.Data.i16[-2]:=(length(h1) div 4)+1;
  X.a:=X.a mod CurveParams.p;
  X.b.Data.i32[length(h1) div 4]:=0;
  for i:=0 to length(h1)-1 do X.b.Data.i8[i]:=h2[i];
  X.b.Data.i16[-1]:=0;
  X.b.Data.i16[-2]:=(length(h1) div 4)+1;
  X.b:=X.b mod CurveParams.p;
  CurrentXp:=0;
  CurrentYp:=0;
  SetPointValue(X);
  if infinity then j:=j+1;
  until not Infinity;
tmp.SetCurveParams(Self.CurveParams);
_Mul_Fp_Fp2Point(CurveParams.Htw,Self,tmp,csJacobian);
Self:=tmp;
end;

{*******************************************************************************}
procedure Fp2Point.SetAsTorsionFromString(s: String);
begin
Self.SetAsTorsionFromHash(SHA256StringHash(s));
end;

{*******************************************************************************}
procedure Fp2Point.SetAsRandomTorsionPoint;
var tmp:Fp2Point;
begin
tmp.SetCurveParams(CurveParams);
Infinity:=false;
Randomize;
repeat
  GetRandomLIntLowerThan(X.a,CurveParams.P);
  GetRandomLIntLowerThan(X.b,CurveParams.P);
  SetPointValue(X);
  if not Infinity then begin
                       if CurveParams.Family=cfBN then _Mul_Fp_Fp2Point(CurveParams.Htw,Self,tmp,csJacobian)
                       else _Fast_Mul_Htw_Fp2Point(CurveParams.u,Self,tmp);
                       Self:=tmp;
                       end;
  until (not Infinity);
  _Mul_Fp_Fp2Point(CurveParams.Rtw,Self,tmp,csAffine);
if not tmp.Infinity then showmessage('problem!');
end;
{*******************************************************************************}
procedure Fp2Point.SetPairingPointCoordinates(PtX, PtY: LInt);
begin
CurrentXp:=PtX;
CurrentYp:=PtY;
end;
{*******************************************************************************}
procedure Fp2Point.SetPointValue(Value: Fp2Int);
begin
_Evaluate_Fp2_Point(Value,Self);
end;

{*******************************************************************************}
procedure Fp2Point.SetToDefaultGenerator;
begin
Self.X:=CurveParams.TwistBasePointX;
Self.Y:=CurveParams.TwistBasePointy;
Self.Z.a:=1;
Self.Z.b:=0;
Self:=CurveParams.Htw*Self;
end;

{*******************************************************************************}
function Fp2Point.ToHexString: String;
begin
if Infinity then Result:='Infinity' else Result:='('+X.ToHexString+','+Y.ToHexString+')';
end;
{*******************************************************************************}
function Fp2Point.ToDecimalString: String;
begin
if Infinity then Result:='Infinity' else Result:='('+X.ToDecimalString+','+Y.ToDecimalString+')';
end;
{*******************************************************************************}
class operator Fp2Point.Add(Left, Right: Fp2Point): Fp2Point;
begin
_Add_Affine_Fp2_Point(Left,Right,Result);
end;
{*******************************************************************************}
class operator Fp2Point.Multiply(Left: LInt; Right{Q}: Fp2Point): Fp2Point;
begin
_Mul_Fp_Fp2Point(Left,Right,Result);
end;
{*******************************************************************************}
class operator Fp2Point.Subtract(Left, Right: Fp2Point): Fp2Point;
begin
_Sub_Fp2_Point(left,right,result);
end;
{*******************************************************************************}
class operator Fp2Point.Negative(const Value: Fp2Point): Fp2Point;
begin
_Neg_Fp2_Point(Value,Result);
end;
{*******************************************************************************}
class operator Fp2Point.Equal(const Left, Right: Fp2Point): Boolean;
begin
Result:=_Are_Equals_Fp2Points(Left,Right);
end;
{*******************************************************************************}
class operator Fp2Point.NotEqual(const Left, Right: Fp2Point): Boolean;
begin
Result:=not _Are_Equals_Fp2Points(Left,Right);
end;

end.

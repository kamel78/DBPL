unit ECCFp;

interface

uses LargeIntegers,Fp2Arithmetic,Fp3Arithmetic,Fp12Arithmetic,Fp6Arithmetic,HashFunctions, GeneralTypes, System.SysUtils,VCL.dialogs;

{*******************************************************************************************
   Développer par Faraoun Kamel Mohamed
   Université Djilali Liabes -Sidi Bel Abbes - Algérie
   kamel_mh@yahoo.fr

   Arithmetic computation over E(Fp): Elliptuc Curves Over Fp.

********************************************************************************************}

Type
   FpPoint=record     //  Definition of a Point on EC Over Fp
      private
      CurveParams: PtrCurveParams;
      public
      X,Y,Z:LInt;
      Lambda,C:LInt;
      LineAtQPFp12:Fp12Int;           // For Tate/Eta Pairings BN Curves
      LineAtQPFp6:Fp6Int;           // For Tate/Eta Pairings MNT Curves
      ComputeLineAtQPFp12:boolean;    // For Opt-Ate/R-ate  BN Curves
      ComputeLineAtQPFp6:boolean;    // For Opt-Ate/R-ate  MNT Curves
      CurrentXQBN,CurrentYQBN:Fp2Int;
      CurrentXQMNT,CurrentYQMNT:Fp3Int;
      Infinity:boolean;
      procedure SetCurveParams(Value:PtrCurveParams);
      procedure SetPointValue(Value:LInt);
      function toDecimalString:String;
      function toHexString:String;
      procedure SetAsRandomPoint;
      procedure SetAsRandomTorsionPoint;
      procedure SetToDefaultGenerator;
      procedure SetAsTorsionFromHash(hash:TBytes);
      procedure SetAsTorsionFromString(s:String);
      function IsOnTheCurve:boolean;
      Procedure JacobianToAffine;
      procedure SetBNPairingPointCoordinates(QtX, QtY: Fp2Int);
      procedure SetMNTPairingPointCoordinates(QtX, QtY: Fp3Int);
      function CompressToArray:TBytes;
      Procedure DeCompressFromArray(a:TBytes);
      class operator Add(const Left, Right: FpPoint): FpPoint;
      class operator Multiply(Left: LInt; const Right: FpPoint): FpPoint;
      class operator Equal(const Left,Right: FpPoint): Boolean;inline;
      class operator NotEqual(const Left, Right: FpPoint): Boolean; inline;
      class operator Negative(const Value: FpPoint): FpPoint;
      end;

  procedure _Add_Jacobian_Fp_Point(left,right:FpPoint;var Result:FpPoint;ComputeLambda:boolean=false);
  procedure _Add_Projective_Fp_Point(Left,Right:FpPoint;var Result:FpPoint;ComputeLambda:boolean=false);
  procedure _Add_Affine_Fp_Point(Left,Right:FpPoint;var Result:FpPoint;ComputeLambda:boolean=false);
  procedure _Double_Jacobian_Fp_Point(const value:FpPoint;var Result:FpPoint;ComputeLambda:boolean=false);
  procedure _Double_Affine_Fp_Point(const Value:FpPoint;var Result:FpPoint;ComputeLambda:boolean=false);
  procedure _Double_Projective_Fp_Point(const Value:FpPoint;var Result:FpPoint;ComputeLambda:boolean=false);
  procedure _Jacobian_To_Affine_FpPoint(Value:FpPoint;var Result:FpPoint);
  procedure _Projective_To_Affine_FpPoint(Value:FpPoint;var Result:FpPoint);
  procedure _Mul_Fp_FpPoint(const Left:LInt;const Right:FpPoint;var Result:FpPoint);
  procedure _Mul_NAF_Fp_FpPoint(const Left:LInt;const Right:FpPoint;var Result:FpPoint);
  procedure _Neg_Fp_Point(Value:FpPoint;var Result:FpPoint);
  procedure _Sub_Fp_Point(Right,Left:FpPoint;var Result:FpPoint);
  function _Are_Equals_FpPoints(left,right:FpPoint):boolean;
  function _Is_OnCurve_FpPoint(value:FpPoint):boolean;


implementation


{*******************************************************************************}
            //      Procedures for Elliptic Curves Arithmetic on FP
{*******************************************************************************}

        {**********   Add two Fp Points using Jacobian coordinates *****************}
        // Left.Z always 1
procedure _Add_Jacobian_Fp_Point(Left,Right:FpPoint;var Result:FpPoint;ComputeLambda:boolean=false);
var t: array[0..10] of LInt;
    Zr2,Zr3:LInt;
begin
Result.CurveParams:=Left.CurveParams;
if ComputeLambda then begin
                      Result.Lambda.Data.i32[-1]:=0;
                      Result.C.Data.i32[-1]:=0;
                      end;
if (Left.Infinity)or(Left.Z.Data.i16[-2]=0) then Result:=Right
else if (Right.Infinity)or(Right.Z.Data.i16[-2]=0) then Result:=Left
else begin
     if _Compare_LInt(Left.X,Right.X)=0 then
                  begin
                  if (_Compare_LInt(Left.Y,Right.Y)=0)and(_Compare_LInt(Left.Z,Right.Z)=0) then _Double_Jacobian_Fp_Point(Left,Result)
                  else begin
                       Result.Infinity:=true;
                       Result.Lambda.Data.i32[-1]:=0;
                      _HCopy_LInt(LEft.X,Result.C);
                      _HCopy_LInt(Left.x,Result.LineAtQPFp12.a.a.a);
                       Result.LineAtQPFp12.a.a.b:=0;
                       Result.LineAtQPFp12.b.a.a:=0;
                       Result.LineAtQPFp12.b.a.b:=0;
                       Result.LineAtQPFp12.a.b:=Right.CurrentXQBN;
                       Result.LineAtQPFp12.a.c.a:=0;
                       Result.LineAtQPFp12.b.b.a:=0;
                       Result.LineAtQPFp12.b.b.b:=0;
                       Result.LineAtQPFp12.a.c.b:=0;
                       Result.LineAtQPFp12.b.c.a:=0;
                       Result.LineAtQPFp12.b.c.b:=0;
                       end;
                  end
     else begin
          Result.Infinity:=false;
          _Sqr_LInt(Right.Z,Zr2);
          _Mod_LInt(Zr2,Left.CurveParams.P,Zr2);
          if ComputeLambda then begin
                                _Mul_LInt(Zr2,Right.Z,t[0]);
                                _Mul_LInt(Left.Y,t[0],t[1]);
                                _Mod_LInt(t[1],Left.CurveParams.P,t[1]);
                                _Sub_LInt(Right.Y,t[1],t[1]);
                                _Mul_LInt(Zr2,Left.X,t[0]);
                                _Mod_LInt(t[0],Left.CurveParams.P,t[0]);
                                _Sub_LInt(Right.X,t[0],t[0]);
                                _Mul_LInt(t[0],Right.Z,t[2]);
                                _Inv_Mod_LInt(t[2],Left.CurveParams.p,t[0]);
                                _Mul_LInt(t[0],t[1],Result.Lambda);
                                _Mod_LInt(Result.Lambda,Left.CurveParams.P,Result.Lambda);
                                _Mul_LInt(Result.Lambda,Left.X,t[0]);
                                _Sub_LInt(Left.Y,t[0],Result.C);
                                _Mod_LInt(Result.C,Left.CurveParams.P,Result.C);
                                end;
          if Result.ComputeLineAtQPFp12 then begin
                                           _Mul_LInt(Zr2,Right.Z,Zr3);
                                           _Mul_LInt(Left.y,Zr3,t[2]);
                                           _Sub_LInt(t[2],Right.Y,t[2]);
                                           _Mod_LInt(t[2],left.CurveParams.P,t[2]);
                                           _Mul_LInt(Left.X,Zr3,t[3]);
                                           _Mul_LInt(Right.X,Right.Z,t[6]);
                                           _Sub_LInt(t[3],t[6],t[3]);
                                           _Mod_LInt(t[3],Left.CurveParams.P,t[3]);
                                           _Mul_LInt(t[3],Right.Y,t[4]);
                                           _Mul_Lint(t[2],Right.x,t[5]);
                                           _Mul_LInt(t[5],Right.Z,t[6]);
                                           _Sub_LInt(t[4],t[6],t[6]);
                                           _Mod_LInt(t[6],Right.CurveParams.P,t[6]);
                                           Result.LineAtQPFp12.a.a.a:=t[6];
                                           Result.LineAtQPFp12.a.a.b:=0;
                                           Result.LineAtQPFp12.b.a.a:=0;
                                           Result.LineAtQPFp12.b.a.b:=0;
                                           _Mul_LInt(Zr3,t[2],t[5]);
                                           _Mod_LInt(t[5],Right.CurveParams.P,t[5]);
                                           _Mul_FP_FP2(t[5],Right.CurrentXQBN,Result.LineAtQPFp12.a.b);
                                           _Mul_LInt(t[3],Zr3,t[5]);
                                           _Mod_LInt(t[5],Right.CurveParams.P,t[5]);
                                           _Sub_LInt(Right.CurveParams.P,t[5],t[5]);
                                           _Mul_FP_FP2(t[5],Right.CurrentYQBN,Result.LineAtQPFp12.b.b);
                                           Result.LineAtQPFp12.a.c.a:=0;
                                           Result.LineAtQPFp12.a.c.b:=0;
                                           Result.LineAtQPFp12.b.c.a:=0;
                                           Result.LineAtQPFp12.b.c.b:=0;
                                           end;
          if Result.ComputeLineAtQPFp6 then begin

                                            end;

          _Mul_LInt(Left.X,Zr2,t[0]);
          _Mod_LInt(t[0],Left.CurveParams.P,t[0]);
          _Sqr_LInt(Left.Y,t[2]);
          _Add_LInt(Left.Y,Right.Z,t[1]);
          _Sqr_LInt(t[1],t[3]);
          _Sub_LInt(t[3],t[2],t[2]);
          _Sub_LInt(t[2],Zr2,t[2]);
          _Mul_LInt(t[2],Zr2,t[3]);
          _Mod_LInt(t[3],Left.CurveParams.P,t[1]);
          _Sub_LInt(t[0],Right.X,t[2]);
          _Mod_LInt(t[2],Left.CurveParams.P,t[2]);
          _Sqr_LInt(t[2],t[3]);
          _Mod_LInt(t[3],Left.CurveParams.P,t[3]);
          _HCopy_LInt(t[3],t[4]);
          _Shl_LInt(t[4],2);
          _Mod_LInt(t[4],Left.CurveParams.P,t[4]);
          _Mul_LInt(t[4],t[2],t[5]);
          _Mod_LInt(t[5],Left.CurveParams.P,t[5]);
          _HCopy_LInt(Right.Y,t[6]);
          _Shl_LInt(t[6],1);
          _Sub_LInt(t[1],t[6],t[6]);
          _Mod_LInt(t[6],Left.CurveParams.P,t[6]);
          _Mul_LInt(t[6],Left.X,t[9]);
          _Mod_LInt(t[9],Left.CurveParams.P,t[9]);
          _Mul_LInt(Right.X,t[4],t[7]);
          _Mod_LInt(t[7],Left.CurveParams.P,t[7]);
          _Shl_LInt(t[7],1);
          _Sqr_LInt(t[6],Result.X);
          _Sub_LInt(Result.X,t[5],Result.X);
          _Sub_LInt(Result.X,t[7],Result.X);
          _Mod_LInt(Result.X,Left.CurveParams.P,Result.X);
          _Shr_LInt(t[7],1);
          _Add_LInt(Right.Z,t[2],Result.Z);
          _Sqr_LInt(Result.Z,t[1]);
          _Sub_LInt(t[1],Zr2,Result.Z);
          _Sub_LInt(Result.Z,t[3],Result.Z);
          _Mod_LInt(Result.Z,Left.CurveParams.P,Result.Z);
          _Sub_LInt(t[7],Result.X,t[1]);
          _Mul_LInt(t[1],t[6],t[8]);
          _Mod_LInt(t[8],Left.CurveParams.P,t[8]);
          _Mul_LInt(Right.Y,t[5],t[0]);
          _Shl_LInt(t[0],1);
          _Mod_LInt(t[0],Left.CurveParams.P,t[0]);
          _Sub_LInt(t[8],t[0],Result.Y);
          _Mod_LInt(Result.Y,Left.CurveParams.P,Result.Y);
          end;
     end;
end;
                {**********   Add two Fp Points using Jacobian coordinates *****************}
                // sould have Right.Z=1 as input every time
                //http://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#addition-madd-1998-cmo
procedure _Add_Projective_Fp_Point(Left,Right:FpPoint;var Result:FpPoint;ComputeLambda:boolean=false);
var t: array[0..9] of LInt;
    t1:Fp3Int;
    tmp:Lint;
begin
Result.CurveParams:=Left.CurveParams;
if ComputeLambda then begin
                      Result.Lambda.Data.i32[-1]:=0;
                      Result.C.Data.i32[-1]:=0;
                      end;
if (Left.Infinity)or(Left.Z.Data.i16[-2]=0) then Result:=Right
else if (Right.Infinity)or(Right.Z.Data.i16[-2]=0) then Result:=Left
else begin
     if _Compare_LInt(Left.X,Right.X)=0 then
                  begin
                  if (_Compare_LInt(Left.Y,Right.Y)=0)and(_Compare_LInt(Left.Z,Right.Z)=0) then _Double_Jacobian_Fp_Point(Left,Result)
                  else begin
                       Result.Infinity:=true;
                       Result.Lambda.Data.i32[-1]:=0;
                      _HCopy_LInt(LEft.X,Result.C);
                      _HCopy_LInt(Left.x,Result.LineAtQPFp12.a.a.a);
                      Result.LineAtQPFp12.a.a.b:=0;
                      Result.LineAtQPFp12.b.a.a:=0; ///w1
                      Result.LineAtQPFp12.b.a.b:=0;
                      Result.LineAtQPFp12.a.b:=Right.CurrentXQBN;
                      Result.LineAtQPFp12.a.c.a:=0;
                      Result.LineAtQPFp12.b.b.a:=0;
                      Result.LineAtQPFp12.b.b.b:=0;
                      Result.LineAtQPFp12.a.c.b:=0;
                      Result.LineAtQPFp12.b.c.a:=0;
                      Result.LineAtQPFp12.b.c.b:=0;
                      if Left.ComputeLineAtQPFp6 then Result.LineAtQPFp6.SetToOne;
                       end;
                  end
     else begin
          Result.Infinity:=false;
          if ComputeLambda then begin     ///SuperSingular Curves
                                _Mul_LInt(Left.Z,Right.Y,t[0]);
                                _Mul_LInt(Left.Z,Right.X,t[1]);
                                _Sub_LInt(t[0],Left.Y,t[3]);
                                _Sub_LInt(t[1],Left.X,t[4]);
                                _Inv_Mod_LInt(t[4],Left.CurveParams.P,t[5]);
                                _Mul_LInt(t[3],t[5],t[6]);
                                _Mod_LInt(t[6],Left.CurveParams.P,Result.Lambda);
                                _Mul_LInt(Result.Lambda,Left.X,t[0]);
                                _Sub_LInt(Left.Y,t[0],t[0]);
                                _Inv_Mod_LInt(Left.Z,Left.CurveParams.P,t[8]);
                                _Mul_LInt(t[8],t[0],Result.C);
                                _Mod_LInt(Result.C,Left.CurveParams.P,Result.C);
                                end;
          _Mul_LInt(Right.Y,Left.Z,t[0]);
          _Sub_LInt(t[0],Left.Y,t[0]);  // u      //*/*/for mnt
          _Sqr_LInt(t[0],t[1]); //u²
          _Mul_LInt(Right.X,Left.Z,t[3]);
          _Sub_LInt(t[3],Left.X,t[3]); //v    //*/*/for mnt
          _Sqr_LInt(t[3],t[4]);  //v²
          _Mul_LInt(t[4],t[3],t[5]);  //v^3
          _Mod_LInt(t[5],Left.CurveParams.P,t[5]);
          _Mul_LInt(t[4],Left.X,t[6]);  //R
          _Mul_LInt(t[1],Left.Z,t[7]);
          _Sub_LInt(t[7],t[5],t[7]);
          _Sub_LInt(t[7],t[6],t[7]);
          _Mod_LInt(t[6],Left.CurveParams.P,t[6]);
          _Mod_LInt(t[7],Left.CurveParams.P,t[7]);
          _Sub_LInt(t[7],t[6],t[7]); //A

      if Left.ComputeLineAtQPFp6 then begin       /// MNT Curves
                                        t1.SetFieldParams(Left.CurveParams.FieldParam);
                                        _Mul_Fp3_FP(Left.CurrentXQMNT,Left.z,t1);
                                        _Sub_Fp3_FP(t1,Left.X,t1);
                                        _Mul_Fp3_FP(t1,t[0],t1);
                                        _Mul_LInt(Left.Y,t[3],tmp);
                                        _Add_Fp3_FP(t1,tmp,t1);
                                       _Neg_Fp3(t1,t1);
                                       _Mul_LInt(Left.Z,t[3],t[1]);
                                       Result.LineAtQPFp6.a.a:=t1.a;
                                       Result.LineAtQPFp6.a.b:=Result.CurrentYQMNT.b*t[1];
                                       Result.LineAtQPFp6.b.a:=Result.CurrentYQMNT.a*t[1];
                                       Result.LineAtQPFp6.b.b:=t1.c;
                                       Result.LineAtQPFp6.c.a:=t1.b;
                                       Result.LineAtQPFp6.c.b:=Result.CurrentYQMNT.c*t[1];
                                         end;
          if Result.ComputeLineAtQPFp12 then begin
                                             _Mul_LInt(t[0],Right.X,t[8]);
                                             _Mul_LInt(t[3],Right.Y,t[9]);
                                             _Sub_LInt(t[8],t[9],Result.LineAtQPFp12.a.a.a);
                                             _Mod_LInt(Result.LineAtQPFp12.a.a.a,Right.CurveParams.P,Result.LineAtQPFp12.a.a.a);
                                             Result.LineAtQPFp12.a.a.b:=0;
                                             _Sub_LInt(Right.CurveParams.P,t[0],t[8]);
                                             _Mul_FP_FP2(t[8],Result.CurrentXQBN,Result.LineAtQPFp12.a.b);
                                             _Mul_FP_FP2(t[3],Result.CurrentYQBN,Result.LineAtQPFp12.b.b);
                                             Result.LineAtQPFp12.b.a.a:=0; ///w1
                                             Result.LineAtQPFp12.b.a.b:=0;
                                             Result.LineAtQPFp12.a.c.a:=0;
                                             Result.LineAtQPFp12.a.c.b:=0;
                                             Result.LineAtQPFp12.b.c.a:=0;
                                             Result.LineAtQPFp12.b.c.b:=0;
                                             end;
          _Mul_LInt(t[7],t[3],Result.X);
          _Sub_LInt(t[6],t[7],t[6]);
          _Mul_LInt(t[6],t[0],Result.Y);
          _Mul_LInt(t[5],Left.Y,t[8]);
          _Sub_LInt(Result.Y,t[8],Result.Y);
          _Mul_LInt(Left.Z,t[5],t[0]);
          _Mod_LInt(Result.X,Left.CurveParams.P,Result.X);
          _Mod_LInt(Result.Y,Left.CurveParams.P,Result.Y);
          _Mod_LInt(t[0],Left.CurveParams.P,Result.Z);
          end;
     end;
end;

            {**********   Add two Fp Points using Affine coordinates *****************}
procedure _Add_Affine_Fp_Point(Left,Right:FpPoint;var Result:FpPoint;ComputeLambda:boolean=false);
var t: array[0..9] of LInt;
    t1:Fp3Int;
begin
Result.CurveParams:=Left.CurveParams;
if ComputeLambda then begin
                      Result.Lambda.Data.i32[-1]:=0;
                      Result.C.Data.i32[-1]:=0;
                      end;
if (Left.Infinity) then Result:=Right
else if (Right.Infinity) then Result:=Left
else begin
     if _Compare_LInt(Left.X,Right.X)=0 then
                  begin
                  if (_Compare_LInt(Left.Y,Right.Y)=0) then _Double_Affine_Fp_Point(Left,Result)
                  else begin
                       Result.Infinity:=true;
                       Result.Lambda.Data.i32[-1]:=0;
                       _HCopy_LInt(LEft.X,Result.C);
                       _HCopy_LInt(Left.x,Result.LineAtQPFp12.a.a.a);
                       Result.LineAtQPFp12.a.a.b:=0;
                       Result.LineAtQPFp12.b.a.a:=0;
                       Result.LineAtQPFp12.b.a.b:=0;
                       Result.LineAtQPFp12.a.b:=Right.CurrentXQBN;
                       Result.LineAtQPFp12.a.c.a:=0;
                       Result.LineAtQPFp12.b.b.a:=0;
                       Result.LineAtQPFp12.b.b.b:=0;
                       Result.LineAtQPFp12.a.c.b:=0;
                       Result.LineAtQPFp12.b.c.a:=0;
                       Result.LineAtQPFp12.b.c.b:=0;
                       if Result.ComputeLineAtQPFp6 then Result.LineAtQPFp6.SetToOne;
                       end;
                  end
     else begin
          Result.Infinity:=false;
          if Result.ComputeLineAtQPFp12 then begin         /// BN Curves Tate,Eta
                                           _Sub_LInt(Left.Y,Right.Y,t[2]);
                                           _Mod_LInt(t[2],left.CurveParams.P,t[2]);
                                           _Sub_LInt(Left.X,Right.X,t[3]);
                                           _Mod_LInt(t[3],Left.CurveParams.P,t[3]);
                                           _Mul_LInt(t[3],Right.Y,t[4]);
                                           _Mul_Lint(t[2],Right.x,t[5]);
                                           _Sub_LInt(t[4],t[5],t[5]);
                                           _Mod_LInt(t[5],Right.CurveParams.P,t[5]);
                                           Result.LineAtQPFp12.a.a.a:=t[5];
                                           Result.LineAtQPFp12.a.a.b:=0;
                                           Result.LineAtQPFp12.b.a.a:=0;
                                           Result.LineAtQPFp12.b.a.b:=0;
                                           _Mul_FP_FP2(t[2],Right.CurrentXQBN,Result.LineAtQPFp12.a.b);
                                           _Sub_LInt(Right.CurveParams.P,t[3],t[5]);
                                           _Mul_FP_FP2(t[5],Right.CurrentYQBN,Result.LineAtQPFp12.b.b);
                                           Result.LineAtQPFp12.a.c.a:=0;
                                           Result.LineAtQPFp12.a.c.b:=0;
                                           Result.LineAtQPFp12.b.c.a:=0;
                                           Result.LineAtQPFp12.b.c.b:=0;
                                           end;
          _Sub_LInt(Right.X,Left.X,t[0]);
          if _IsNeg(t[0]) then _Add_LInt(t[0],Left.CurveParams.P,t[0]);
          _Sub_LInt(Right.Y,Left.Y,t[1]);
          if _IsNeg(t[1]) then _Add_LInt(t[1],Left.CurveParams.P,t[1]);
          _Inv_Mod_LInt(t[0],Left.CurveParams.P,t[0]);
          _Mul_LInt(t[1],t[0],t[2]);     //// Lamda
          _Sqr_LInt(t[2],t[3]);
          _Sub_LInt(t[3],Left.X,t[3]);
          _Sub_LInt(t[3],Right.X,t[3]);
          if Result.ComputeLineAtQPFp6 then begin   //// MNT Curves
                                            t1.SetFieldParams(Left.CurveParams.FieldParam);
                                            _Sub_Fp3_FP(Right.CurrentXQMNT,Right.X,t1);
                                            _Mul_FP_Fp3(t[2],t1,t1);
                                            _Add_FP_Fp3(Right.Y,t1,t1);
                                            _Neg_Fp3(t1,t1);
                                            Result.LineAtQPFp6.a.a:=t1.a;
                                            Result.LineAtQPFp6.a.b:=Right.CurrentYQMNT.b;
                                            Result.LineAtQPFp6.b.a:=Right.CurrentYQMNT.a;
                                            Result.LineAtQPFp6.b.b:=t1.c;
                                            Result.LineAtQPFp6.c.a:=t1.b;
                                            Result.LineAtQPFp6.c.b:=Right.CurrentYQMNT.c;
                                            end;
          _Mod_LInt(t[3],Left.CurveParams.P,Result.X);
          _Sub_LInt(Left.X,Result.X,t[4]);
          if _IsNeg(t[0]) then _Add_LInt(t[4],Left.CurveParams.P,t[4]);
          _Mul_LInt(t[2],t[4],Result.Y);
          _Sub_LInt(Result.Y,Left.Y,Result.Y);
          _Mod_LInt(Result.Y,Left.CurveParams.P,Result.Y);
          if ComputeLambda then begin                // SuperSingular Curves
                                _HCopy_LInt(t[2],Result.Lambda);
                                _Mul_LInt(t[2],Left.X,t[5]);
                                _Sub_LInt(Left.Y,t[5],Result.C);
                                _Mod_LInt(Result.C,Left.CurveParams.P,Result.C);
                                end;
          end;
     end;
end;

        {**********   Double an Fp Point using Jacobian coordinates  *****************}
procedure _Double_Jacobian_Fp_Point(const Value:FpPoint;var Result:FpPoint;ComputeLambda:boolean=false);
var t: array[0..6] of LInt;
    ZQ2:LInt;
begin
Result.CurveParams:=Value.CurveParams;
if ComputeLambda then begin
                      Result.Lambda.Data.i32[-1]:=0;
                      Result.C.Data.i32[-1]:=0;
                      end;
if(Value.Infinity)or(Value.Y.Data.i16[-2]=0)or(Value.Z.Data.i16[-2]=0) then Result.Infinity:=true
else begin
     Result.Infinity:=false;
     _Sqr_LInt(Value.Z,ZQ2);
     _Mod_LInt(ZQ2,Value.CurveParams.P,ZQ2);
     _Sqr_LInt(Value.X,t[0]);
     _Mod_LInt(t[0],Value.CurveParams.P,t[0]);
     if ComputeLambda then begin
                           if Value.CurveParams.A.Data.i16[-2]<>0 then begin
                                                                    _Sqr_LInt(ZQ2,t[1]);
                                                                    _Mul_LInt(t[1],Value.CurveParams.A,t[2]);
                                                                    end;
                           _HCopy_LInt(t[2],t[1]); //z^4
                           _HCopy_LInt(t[0],t[2]);  //x^2
                           _Shl_LInt(t[2],1);
                           _Add_LInt(t[0],t[2],t[2]);
                           _Add_LInt(t[2],t[1],t[2]);
                           _Mod_LInt(t[2],Value.CurveParams.P,t[2]);  //3x^2
                           _Inv_Mod_LInt(Value.Z,Value.CurveParams.P,t[4]);
                           _Inv_Mod_LInt(Value.Y,Value.CurveParams.P,t[5]);
                           _Mul_LInt(t[4],t[5],t[6]);
                           _Mul_LInt(t[6],Value.CurveParams.FieldParam.inv_2,t[3]);
                           _Mul_LInt(t[2],t[3],Result.Lambda);
                           _Mod_LInt(Result.Lambda,Value.CurveParams.P,Result.Lambda);
                           _Mul_LInt(Result.Lambda,Value.X,t[1]);
                           _Mul_LInt(t[1],Value.Z,t[2]);
                           _Sub_LInt(Value.Y,t[2],Result.C);
                           _Sqr_LInt(t[4],t[5]);
                           _Mul_LInt(t[5],t[4],t[2]);
                           _Mul_LInt(t[2],Result.C,t[4]);
                           _Mod_LInt(t[4],Value.CurveParams.P,Result.C);
                           end;
     if Value.ComputeLineAtQPFp12 then begin
                                       ///  to be used for Eta pairing on BN curves
                                       _HCopy_LInt(t[0],t[2]);
                                       _Shl_LInt(t[2],1);
                                       _Add_LInt(t[0],t[2],t[2]);
                                       _Mod_LInt(t[2],Value.CurveParams.P,t[2]);  // n=3x^2
                                       _Mul_LInt(Value.y,Value.Z,t[3]);
                                       _Shl_LInt(t[3],1);
                                       _Mod_LInt(t[3],Value.CurveParams.P,t[3]);  // d=2y*z
                                       _Mul_LInt(t[3],Value.Y,t[4]);
                                       _Mul_Lint(t[2],Value.x,t[5]);
                                       _Mul_LInt(t[5],Value.Z,t[6]);
                                       _Sub_LInt(t[4],t[6],t[6]);
                                       _Mod_LInt(t[6],value.CurveParams.P,t[6]);
                                       Result.LineAtQPFp12.a.a.a:=t[6];
                                       Result.LineAtQPFp12.a.a.b:=0;
                                       Result.LineAtQPFp12.b.a.a:=0; ///w1
                                       Result.LineAtQPFp12.b.a.b:=0;
                                       _Mul_LInt(Zq2,Value.Z,t[4]);   //// Z^3
                                       _Mul_LInt(t[4],t[2],t[5]);
                                       _Mod_LInt(t[5],Value.CurveParams.P,t[5]);
                                       _Mul_FP_FP2(t[5],value.CurrentXQBN,Result.LineAtQPFp12.a.b);
                                       _Mul_LInt(t[3],t[4],t[5]);
                                       _Mod_LInt(t[5],Value.CurveParams.P,t[5]);
                                       _Sub_LInt(Value.CurveParams.P,t[5],t[5]);
                                       _Mul_FP_FP2(t[5],value.CurrentYQBN,Result.LineAtQPFp12.b.b);
                                       Result.LineAtQPFp12.a.c.a:=0;
                                       Result.LineAtQPFp12.a.c.b:=0;
                                       Result.LineAtQPFp12.b.c.a:=0;
                                       Result.LineAtQPFp12.b.c.b:=0;
                                       end;
     _Sqr_LInt(Value.Y,t[1]);
     _Mod_LInt(t[1],Value.CurveParams.P,t[1]);
     _Sqr_LInt(t[1],t[2]);
     _Mod_LInt(t[2],Value.CurveParams.P,t[2]);
     _Add_LInt(t[1],Value.X,t[3]);
     _Sqr_LInt(t[3],t[4]);
     _Sub_LInt(t[4],t[0],t[3]);
     _Sub_LInt(t[3],t[2],t[3]);
     _Shl_LInt(t[3],1);
     _Mod_LInt(t[3],Value.CurveParams.P,t[3]);
     _HCopy_LInt(t[0],t[4]);
     _Shl_LInt(t[4],1);
     _Add_LInt(t[4],t[0],t[4]);
     if Value.CurveParams.A.Data.i16[-2]<>0 then begin
                                        _Sqr_LInt(ZQ2,t[6]);
                                        _Mul_LInt(Value.CurveParams.A,t[6],t[5]);
                                        _Add_LInt(t[4],t[5],t[4]);
                                        end;
     _Mod_LInt(t[4],Value.CurveParams.P,t[4]);
     _Add_LInt(Value.X,t[4],t[6]);
     _Mod_LInt(t[6],Value.CurveParams.P,t[6]);
     _Sqr_LInt(t[4],t[5]);
     _Mod_LInt(t[5],Value.CurveParams.P,t[5]);
     _Shl_LInt(t[3],1);
     _Sub_LInt(t[5],t[3],t[5]);
     _Mod_LInt(t[5],Value.CurveParams.P,Result.X);
     _Shr_LInt(t[3],1);
     _Add_LInt(Value.Y,Value.Z,Result.Z);
     _Sqr_LInt(Result.Z,t[0]);
     _Sub_LInt(t[0],t[1],Result.Z);
     _Sub_LInt(Result.Z,ZQ2,Result.Z);
     _Mod_LInt(Result.Z,Value.CurveParams.P,Result.Z);
     _Sub_LInt(t[3],Result.X,Result.Y);
     _Mul_LInt(Result.Y,t[4],t[0]);
     _Shl_LInt(t[2],3);
     _Sub_LInt(t[0],t[2],Result.Y);
     _Mod_LInt(Result.Y,Value.CurveParams.P,Result.Y);
     end;
end;

        {**********   Double an Fp Point using Projective coordinates (X/Z,Y/Z) *****************}
        //     http://www.hyperelliptic.org/EFD/g1p/auto-shortw-projective.html#doubling-dbl-2007-bl
procedure _Double_Projective_Fp_Point(const Value:FpPoint;var Result:FpPoint;ComputeLambda:boolean=false);
var t: array[0..8] of LInt;
    t1:Fp3Int;
    a1,a2,a3:Lint;
begin
Result.CurveParams:=Value.CurveParams;
if ComputeLambda then begin
                      Result.Lambda.Data.i32[-1]:=0;
                      Result.C.Data.i32[-1]:=0;
                      end;
if Value.ComputeLineAtQPFp12 then Value.LineAtQPFp12.SetToOne;
if(Value.Infinity)or(Value.Y.Data.i16[-2]=0)or(Value.Z.Data.i16[-2]=0) then Result.Infinity:=true
else begin
     Result.Infinity:=false;
     _Sqr_LInt(value.X,t[0]);
     _Add_LInt(t[0],t[0],t[1]);
     _Add_LInt(t[0],t[1],t[1]);
      if ComputeLambda then begin
                           _Sqr_LInt(Value.Z,t[2]);
                           _Mul_LInt(t[2],Value.CurveParams.A,t[3]);
                           _Add_LInt(t[1],t[3],t[3]);
                           _Mul_LInt(Value.Z,Value.Y,t[4]);
                           _Shl_LInt(t[4],1);
                           _Inv_Mod_LInt(t[4],Value.CurveParams.P,t[5]);
                           _Mul_LInt(t[5],t[3],Result.Lambda);
                           _Mod_LInt(Result.Lambda,Value.CurveParams.P,Result.Lambda);
                           _Mul_LInt(Result.Lambda,Value.X,t[6]);
                           _Sub_LInt(Value.Y,t[6],t[7]);
                           _Inv_Mod_LInt(Value.Z,Value.CurveParams.P,t[8]);
                           _Mul_LInt(t[8],t[7],Result.C);
                           _Mod_LInt(Result.C,Value.CurveParams.P,Result.C);
                           end;
      _Sqr_LInt(Value.Z,t[2]);
      _Mul_LInt(t[2],Value.CurveParams.A,t[3]);
      _Add_LInt(t[1],t[3],t[1]);      //3*x2+a*z²
      _Mul_LInt(value.Y,value.Z,t[2]);
      _Shl_LInt(t[2],1);
      if Value.ComputeLineAtQPFp6 then begin              /// MNT Curves
                                        t1.SetFieldParams(Value.CurveParams.FieldParam);
                                        _Mul_Fp3_FP(Value.CurrentXQMNT,Value.Z,t1);
                                        _sub_Fp3_FP(t1,Value.X,t1);
                                        _Mul_Fp_FP3(t[1],t1,t1);
                                        _Mul_LInt(t[2],value.y,t[3]);
                                        _Mul_LInt(t[2],value.z,t[4]);
                                        _add_Fp3_FP(t1,t[3],t1);
                                       _Neg_Fp3(t1,t1);
                                       Result.LineAtQPFp6.a.a:=t1.a;
                                       Result.LineAtQPFp6.a.b:=Value.CurrentYQMNT.b*t[4];
                                       Result.LineAtQPFp6.b.a:=Value.CurrentYQMNT.a*t[4];
                                       Result.LineAtQPFp6.b.b:=t1.c;
                                       Result.LineAtQPFp6.c.a:=t1.b;
                                       Result.LineAtQPFp6.c.b:=Value.CurrentYQMNT.c*t[4];
                                       end;



      _Sqr_LInt(t[2],t[3]);
       if Value.ComputeLineAtQPFp12 then begin
                                          _Sqr_LInt(Value.Z,t[4]);
                                          _Sqr_LInt(Value.Y,t[5]);
                                          _Add_LInt(t[4],t[4],t[6]);
                                          _Add_LInt(t[6],t[4],t[4]);
                                          _Mul_LInt(Value.CurveParams.B,t[4],t[6]);
                                          _Sub_LInt(t[6],t[5],Result.LineAtQPFp12.a.a.a);
                                          Result.LineAtQPFp12.a.a.b:=0;
                                          _Mod_LInt(Result.LineAtQPFp12.a.a.a,Value.CurveParams.P,Result.LineAtQPFp12.a.a.a);
                                          _Mul_FP_FP2(t[1],Result.CurrentXQBN,Result.LineAtQPFp12.a.b);
                                          _Sub_LInt(Value.CurveParams.P,t[2],t[6]);
                                          _Mul_FP_FP2(t[6],Result.CurrentYQBN,Result.LineAtQPFp12.b.b);
                                          Result.LineAtQPFp12.b.a.a:=0;
                                          Result.LineAtQPFp12.b.a.b:=0;
                                          Result.LineAtQPFp12.a.c.a:=0;
                                          Result.LineAtQPFp12.a.c.b:=0;
                                          Result.LineAtQPFp12.b.c.a:=0;
                                          Result.LineAtQPFp12.b.c.b:=0;
                                          end;
      _Mul_LInt(t[3],t[2],Result.Z);
      _Mul_LInt(t[2],Value.Y,t[4]);
      _Sqr_LInt(t[4],t[8]);
      _Add_LInt(Value.X,t[4],t[5]);
      _Sqr_LInt(t[5],t[6]);
      _Sub_LInt(t[6],t[0],t[6]);
      _Sub_LInt(t[6],t[8],t[6]);
      _Sqr_LInt(t[1],t[0]);
      _Sub_LInt(t[0],t[6],t[0]);
      _Sub_LInt(t[0],t[6],t[0]);
      _Mul_LInt(t[2],t[0],t[7]);
      _Mod_LInt(t[7],Value.CurveParams.P,t[7]);
      _HCopy_LInt(t[7],Result.X);
      _Sub_LInt(t[6],t[0],t[6]);
      _Mul_LInt(t[1],t[6],Result.Y);
      _Sub_LInt(Result.Y,t[8],Result.Y);
      _Sub_LInt(Result.Y,t[8],Result.Y);
      _Mod_LInt(Result.Y,Value.CurveParams.P,Result.Y);
      _Mod_LInt(Result.Z,Value.CurveParams.P,Result.Z);
      end;
end;

        {**********   Double an Fp Point using Affine coordinates  *****************}
procedure _Double_Affine_Fp_Point(const Value:FpPoint;var Result:FpPoint;ComputeLambda:boolean=false);
var t: array[0..7] of LInt;
    t1:Fp3Int;
begin
Result.CurveParams:=Value.CurveParams;
if ComputeLambda then begin
                      Result.Lambda.Data.i32[-1]:=0;
                      Result.C.Data.i32[-1]:=0;
                      end;
if(Value.Infinity)or(Value.Y.Data.i16[-2]=0) then Result.Infinity:=true
else begin
     Result.Infinity:=false;
     _Sqr_LInt(Value.X,t[0]);
     _Add_LInt(t[0],t[0],t[1]);
     _Add_LInt(t[0],t[1],t[1]);
     _Add_LInt(t[1],Value.CurveParams.A,t[1]);
     if Value.ComputeLineAtQPFp12 then begin     //BN Curves
                                       ///  to be used for Tate and Eta pairing
                                       _Sqr_LInt(Value.Y,t[3]);
                                       _Shl_LInt(t[3],1);
                                       _Mul_LInt(t[1],Value.X,t[2]);
                                       _Sub_LInt(t[3],t[2],t[3]);
                                       _Mod_LInt(t[3],Value.CurveParams.P,t[3]);
                                       Result.LineAtQPFp12.a.a.a:=t[3];
                                       Result.LineAtQPFp12.a.a.b:=0;
                                       Result.LineAtQPFp12.b.a.a:=0; ///w1
                                       Result.LineAtQPFp12.b.a.b:=0;
                                       _Mul_FP_FP2(t[1],Value.CurrentXQBN,Result.LineAtQPFp12.a.b);
                                       _Add_LInt(Value.Y,Value.Y,t[2]);
                                       _Sub_LInt(Value.CurveParams.P,t[2],t[2]);
                                       _Mul_FP_FP2(t[2],Value.CurrentYQBN,Result.LineAtQPFp12.b.b);
                                       Result.LineAtQPFp12.a.c.a:=0;
                                       Result.LineAtQPFp12.a.c.b:=0;
                                       Result.LineAtQPFp12.b.c.a:=0;
                                       Result.LineAtQPFp12.b.c.b:=0;
                                       end;
     _Add_LInt(Value.Y,Value.Y,t[2]);
     _Inv_Mod_LInt(t[2],Value.CurveParams.P,t[3]);
     _Mul_LInt(t[1],t[3],t[4]);
     _Mod_LInt(t[4],Value.CurveParams.P,t[4]); //Lambda
     _Sqr_LInt(t[4],t[5]);
     _Sub_LInt(t[5],Value.X,t[5]);
     _Sub_LInt(t[5],Value.X,t[5]);
     _Mod_LInt(t[5],Value.CurveParams.P,t[5]);
     _Sub_LInt(Value.X,t[5],t[6]);
     _Mul_LInt(t[4],t[6],t[7]);
     _Sub_LInt(t[7],Value.Y,t[7]);
     if Value.ComputeLineAtQPFp6  then begin     //MNT Curves
                                       t1.SetFieldParams(Value.CurveParams.FieldParam);
                                       _Sub_Fp3_FP(Value.CurrentXQMNT,Value.X,t1);
                                       _Mul_FP_Fp3(t[4],t1,t1);
                                       _Add_FP_Fp3(Value.Y,t1,t1);
                                       _Neg_Fp3(t1,t1);
                                       Result.LineAtQPFp6.a.a:=t1.a;
                                       Result.LineAtQPFp6.a.b:=Value.CurrentYQMNT.b;
                                       Result.LineAtQPFp6.b.a:=Value.CurrentYQMNT.a;
                                       Result.LineAtQPFp6.b.b:=t1.c;
                                       Result.LineAtQPFp6.c.a:=t1.b;
                                       Result.LineAtQPFp6.c.b:=Value.CurrentYQMNT.c;
                                       end;

     if ComputeLambda then begin               // SuperSingular Curves
                           _HCopy_LInt(t[4],Result.Lambda);
                           _Mul_LInt(t[4],Value.X,t[6]);
                           _Sub_LInt(Value.Y,t[6],Result.C);
                           _Mod_LInt(Result.C,Value.CurveParams.P,Result.C);
                           end;
     _Mod_LInt(t[7],Value.CurveParams.P,Result.Y);
     _HCopy_LInt(t[5],Result.X);

     end;
end;
         {**********   Multiply a scalar with an Fp Point *******}
         //Right Should be different than Result
procedure _Mul_Fp_FpPoint(const Left:LInt;const Right:FpPoint;var Result:FpPoint);
var i:integer;
begin
Result.CurveParams:=Right.CurveParams;
if (Right.Infinity)or (Left.Data.i16[-2]=0) then Result.Infinity:=true
else begin
     Result.Infinity:=true;
     if left=1 then Result:=Right
     else if Left=2 then _Double_Jacobian_Fp_Point(Right,Result)
     else
     for i:=Left.BitLength-1 downto 0 do begin
                                         _Double_Jacobian_Fp_Point(Result,Result);
                                         if _Is_BitSet_At(Left,i) then
                                                  _Add_Jacobian_Fp_Point(Right,Result,Result);
                                         end;
     _Jacobian_To_Affine_FpPoint(Result,Result);
     if _IsNeg(Left) then _Neg_Fp_Point(Result,Result);
     end;
end;
       {**********   Multiply a scalar with an Fp Point *******}
       //Right Should be different than Result
procedure _Mul_NAF_Fp_FpPoint(const Left:LInt;const Right:FpPoint;var Result:FpPoint);
var i:integer;
    loop:LIntArrayForm;
    NegLeft:FpPoint;
begin
Result.CurveParams:=Right.CurveParams;
_Neg_Fp_Point(Right,NegLeft);
if (Right.Infinity)or (Left.Data.i16[-2]=0) then Result.Infinity:=true
else begin
     Result.Infinity:=true;
     Loop:=Left.ToNafArray;
     for i:=Length(Loop)-1 downto 0 do begin
                                         _Double_Jacobian_Fp_Point(Result,Result);
                                         if Loop[i]=1 then
                                                _Add_Jacobian_Fp_Point(Right,Result,Result);
                                         if Loop[i]=-1 then
                                                _Add_Jacobian_Fp_Point(NegLeft,Result,Result);
                                         end;
     _Jacobian_To_Affine_FpPoint(Result,Result);
     end;
end;

       {**********   Convert Jacobian to Affine for an Fp Point *****}
procedure _Jacobian_To_Affine_FpPoint(Value:FpPoint;var Result:FpPoint);
var t1,t2:LInt;
begin
if (Value.Infinity) or (Value.Z.Data.i16[-2]=0) then Result.Infinity:=true
else begin
     t1:=Value.z.InversModulo(Value.CurveParams.P);
     t2:=t1.Sqr mod Value.CurveParams.P;
     Result.X:=(Value.X*t2) mod Value.CurveParams.P;
     Result.y:=(Value.y*t1) mod Value.CurveParams.P;
     Result.y:=(Result.y*t2) mod Value.CurveParams.P;
     Result.z:=1;
     end;
end;
       {**********   Convert Projective to Affine for an Fp Point *****}
procedure _Projective_To_Affine_FpPoint(Value:FpPoint;var Result:FpPoint);
var t1:LInt;
begin
if (Value.Infinity) or (Value.Z.Data.i16[-2]=0) then Result.Infinity:=true
else begin
     _Inv_Mod_LInt(Value.Z,Value.CurveParams.P,t1);
     _Mul_LInt(Value.x,t1,Result.x);
     _Mod_LInt(Result.X,Value.CurveParams.P,Result.X);
     _Mul_LInt(Value.Y,t1,Result.Y);
     _Mod_LInt(Result.y,Value.CurveParams.P,Result.y);
     Result.z:=1;
     end;
end;

        {**********   Negation of an Fp Point *************}
procedure _Neg_Fp_Point(Value:FpPoint;var Result:FpPoint);
begin
  Result:=Value;
  _Sub_LInt(Value.CurveParams.P,Result.Y,Result.Y);
end;
        {**********   Substract two Fp Points **************}
procedure _Sub_Fp_Point(Right,Left:FpPoint;var Result:FpPoint);
var tmp:FpPoint;
begin
  _Neg_Fp_Point(Right,tmp);
  _Add_Jacobian_Fp_Point(tmp,Left,Result);
end;
      {**********   Evaluate an Fp Point from X value ******}
procedure _Evaluate_Fp_Point(value:LInt; var Result:FpPoint);
var tmp,tmp1:LInt;
begin
Result.X:=Value;
Result.Z:=1;
Result.Infinity:=false;
_Sqr_LInt(Value,tmp);
_Mod_LInt(tmp,Result.CurveParams.P,tmp);
_Add_LInt(tmp,Result.CurveParams.A,tmp);
_Mul_LInt(tmp,Value,tmp1);
_Add_LInt(tmp1,Result.CurveParams.B,tmp);
_Mod_LInt(tmp,Result.CurveParams.P,tmp);
if tmp.IsASqrMod(Result.CurveParams.P,Result.CurveParams.FieldParam.MontgomeryData)
       then begin
            ModSquareLInt(tmp,Result.CurveParams.P, result.Y,Result.CurveParams.FieldParam.MontgomeryData,false);
             Result.Y:=Result.CurveParams.P-Result.Y;
            end
else Result.Infinity:=true;
end;
       {**********   Test if two Fp Points are equal *******}
function _Are_Equals_FpPoints(left,right:FpPoint):boolean;
begin
Result:=(Left.X=Right.X)and(Left.Y=Right.Y);
end;
   {**********   Test if an Fp Point is on the Curve ***********}
function _Is_OnCurve_FpPoint(value:FpPoint):boolean;
var tmp,tmp1:LInt;
begin
_Sqr_LInt(Value.X,tmp);
_Add_LInt(tmp,Value.CurveParams.A,tmp);
_Mul_LInt(tmp,Value.X,tmp1);
_Add_LInt(tmp1,Value.CurveParams.B,tmp);
_Mod_LInt(tmp,Value.CurveParams.P,tmp);
_Sqr_LInt(Value.Y,tmp1);
_Mod_LInt(tmp1,Value.CurveParams.P,tmp1);
Result:=(tmp.IsASqrMod(Value.CurveParams.P,Value.CurveParams.FieldParam.MontgomeryData)) and (tmp1=tmp);
end;

{*******************************************************************************}
      ///  Definitions of an Fp Point operators and functions
{*******************************************************************************}
procedure FpPoint.SetPointValue(Value:LInt);
begin
_Evaluate_Fp_Point(value,Self);
end;

{*******************************************************************************}
procedure FpPoint.SetToDefaultGenerator;
var tmp:FpPoint;
begin
tmp.SetCurveParams (Self.CurveParams);
tmp.X:=CurveParams.BasePointX;
tmp.Y:=CurveParams.BasePointY;
tmp.Z:=1;
_Mul_Fp_FpPoint(CurveParams.H,tmp,Self);
end;

{*******************************************************************************}
procedure FpPoint.SetCurveParams(Value:PtrCurveParams);
begin
CurveParams:=Value;
Infinity:=false;
ComputeLineAtQPFp12:=false;
Z:=1;
end;
{*******************************************************************************}
function FpPoint.toDecimalString:String;
begin
if Infinity then Result:='Infinity'
else Result:='('+X.ToDecimalString+' , '+y.ToDecimalString+')';
end;
{*******************************************************************************}
function FpPoint.toHexString:String;
begin
if Infinity then Result:='Infinity'
else Result:='('+X.ToHexString+' , '+Y.ToHexString+')';
end;
{*******************************************************************************}
function FpPoint.IsOnTheCurve: boolean;
begin
Result:=_Is_OnCurve_FpPoint(Self);
end;
{*******************************************************************************}
procedure FpPoint.JacobianToAffine;
begin
_Jacobian_To_Affine_FpPoint(Self,Self);
end;
{*******************************************************************************}
procedure FpPoint.SetBNPairingPointCoordinates(QtX, QtY: Fp2Int);
begin
CurrentXQBN:=QtX;
CurrentYQBN:=Qty;
end;
{*******************************************************************************}
procedure FpPoint.SetMNTPairingPointCoordinates(QtX, QtY: Fp3Int);
begin
CurrentXQMNT:=QtX;
CurrentYQMNT:=Qty;
end;

{*******************************************************************************}
procedure FpPoint.SetAsRandomPoint;
begin
Infinity:=False;
Randomize;
Repeat
   GetRandomLIntLowerThan(Self.X,CurveParams.P);
   _Evaluate_Fp_Point(X,Self);
   _Mod_LInt(Y,CurveParams.P,Y);
  Until (not Infinity) and( Y.Data.i16[-2]<>0);
end;
{*******************************************************************************}
procedure FpPoint.SetAsTorsionFromHash(hash: TBytes);
var len,i:integer;
    tmppoint:FpPoint;
begin
len:=length(hash);
X.Data.i32[len div 4]:=0;
for i:=0 to len-1 do X.Data.i8[i]:=hash[i];
X.Data.i16[-1]:=0;
X.Data.i16[-2]:=(len div 4)+1;
Infinity:=false;
X:=X mod CurveParams.P;
repeat
  SetPointValue(X);
  if not Infinity then _Mul_Fp_FpPoint(CurveParams.H,Self,Tmppoint)
  else _Inc_LInt(X,1);
  until not Infinity;
Self:=Tmppoint;
end;

{*******************************************************************************}
procedure FpPoint.SetAsTorsionFromString(s:String);
begin
Self.SetAsTorsionFromHash(SHA256StringHash(s));
end;

{*******************************************************************************}
procedure FpPoint.SetAsRandomTorsionPoint;
var tmp:FpPoint;
begin
Infinity:=false;
tmp.SetCurveParams(Self.CurveParams);
Randomize;
repeat
  GetRandomLIntLowerThan(X,CurveParams.P);
  SetPointValue(X);
  if not Infinity then begin
                       _Mul_Fp_FpPoint(CurveParams.H,Self,tmp);
                       Self:=tmp;
                       end;
until (not Infinity);
{  _Mul_Fp_FpPoint(CurveParams.R,Self,tmp);
if not tmp.Infinity then showmessage('problem! G1');}
end;
{*******************************************************************************}
function FpPoint.CompressToArray:TBytes;
var i:Integer;
begin
if X.BitLength mod 8=0 then Setlength(Result,(X.BitLength div 8))
else Setlength(Result,(X.BitLength div 8)+1);
SetLength(Result,Length(Result)+1);
for i:=0 to Length(Result)-2 do Result[i]:=X.Data.i8[i];
if (Y.Data.i8[0] and 1)=0 then Result[Length(Result)-1]:=0
else Result[Length(Result)-1]:=1;
end;
{*******************************************************************************}
procedure FpPoint.DeCompressFromArray(a:TBytes);
var i:Integer;
    sign:byte;
    tmp:LInt;
begin
Sign:=a[length(a)-1];
X.Data.i16[-1]:=0;
X.Data.i32[((Length(a)-1) div 4)]:=0;
for i:=0 to length(a)-2 do X.Data.i8[i]:=a[i];
X.Data.i16[-2]:=((Length(a)-1) div 4)+1;
tmp:=((X.Sqr*X)+(CurveParams.A*X)+ CurveParams.B) mod CurveParams.P;
if tmp.IsASqrMod(CurveParams.P,CurveParams.FieldParam.MontgomeryData)
       then begin
            ModSquareLInt(tmp,CurveParams.P, Y,CurveParams.FieldParam.MontgomeryData,false);
            if (Y.Data.i8[0] and 1)<>sign then _Sub_LInt(CurveParams.P,Y,Y);
            end
else Infinity:=true;
end;
{*******************************************************************************}
class operator FpPoint.NotEqual(const Left, Right: FpPoint): Boolean;
begin
Result:=not(_Are_Equals_FpPoints(left,right));
end;
{*******************************************************************************}
class operator FpPoint.Negative(const Value: FpPoint): FpPoint;
begin
_Neg_Fp_Point(Value,Result);
end;
{*******************************************************************************}
class operator FpPoint.Add(const Left, Right: FpPoint): FpPoint;
begin
_Add_Jacobian_Fp_Point(left,right,Result);
_Jacobian_To_Affine_FpPoint(Result,Result);
end;
{*******************************************************************************}
class operator FpPoint.Multiply(Left: LInt; const Right: FpPoint): FpPoint;
begin
_Mul_Fp_FpPoint(Left,Right,Result);
end;
{*******************************************************************************}
class operator FpPoint.Equal(const Left,Right: FpPoint): Boolean;
begin
Result:=_Are_Equals_FpPoints(Left,Right)
end;
end.

using Symbolics
using GenericLinearAlgebra

@variables s1 i1 d1 t1 s2 i2 d2 t2 s3 i3 d3 t3 s4 i4 d4 t4 
@variables α λ1 γ1 δ1 σ1 τ1     # group 1 parameters 
@variables λ2 γ2 δ2 σ2 τ2       # group 2 parameters 
@variables λ3 γ3 δ3 σ3 τ3       # group 3 parameters 
@variables λ4 γ4 δ4 σ4 τ4       # group 4 parameters 
@variables c11 c12 c13 c14   c21 c22 c23 c24   c31 c32 c33 c34    c41 c42 c43 c44   # contact Matrix entry

F = [
        -s1 * α * (c11 * i1 + c12 * i2 + c13 * i3 + c14 * i4);
        0;
        0;
        
        -s2 * α * (c21 * i1 + c22 * i2 + c23 * i3 + c24 * i4);
        0;
        0;
        
        -s3 * α * (c31 * i1 + c32 * i2 + c33 * i3 + c34 * i4);
        0;
        0;
        
        -s4 * α * (c41 * i1 + c42 * i2 + c43 * i3 + c44 * i4);
        0;
        0 ]
        
V = [ ((λ1*γ1)/(λ1+γ1)) * i1  + γ1 * i1;
        i1 * γ1 - d1 * (λ1 + δ1); 
        δ1 * d1 - (σ1 + τ1) * t1;

        ((λ2*γ2)/(λ2+γ2)) * i2  + γ2 * i2;
        i2 * γ2 - d2 * (λ2 + δ2); 
        δ2 * d2 - (σ1 + τ1) * t1;

        ((λ3*γ3)/(λ3+γ3)) * i3  + γ3 * i3;
        i3 * γ3 - d3 * (λ3 + δ3); 
        δ3 * d3 - (σ3 + τ3) * t3;

        ((λ4*γ4)/(λ4+γ4)) * i4  + γ4 * i4;
        i4 * γ4 - d4 * (λ4 + δ4); 
        δ4 * d4 - (σ4 + τ4) * t4 ]

X = [ i1, d1, t1, i2, d2, t2, i3, d3, t3, i4, d4, t4 ]


JacF = Symbolics.jacobian(F, X)
JacV = Symbolics.jacobian(V, X)

substitute(JacF, Dict([i1 => 0, d1 => 0, t1 => 0, i2 => 0, d2 => 0, t2 => 0, i3 => 0, d3 => 0, t3 => 0, i4 => 0, d4 => 0, t4 => 0]) )
substitute(JacV, Dict([i1 => 0, d1 => 0, t1 => 0, i2 => 0, d2 => 0, t2 => 0, i3 => 0, d3 => 0, t3 => 0, i4 => 0, d4 => 0, t4 => 0]) )

A = JacF*inv(JacV) 
GenericLinearAlgebra.eigvals(A)
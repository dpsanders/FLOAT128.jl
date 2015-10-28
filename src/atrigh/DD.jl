# |a| <= 0.00275
function atanh_series(a::DD)
    x = a
    x2 = x*x
    x4 = x2*x2
    x6 = x2*x4
    x8 = x4*x4
    x16 = x8*x8
    x24 = x8*x16

    z = dd_recip_int[3]*x2 + dd_recip_int[5]*x4 + dd_recip_int[7]*x6
    z1 = x8*(dd_recip_int[9] + dd_recip_int[11]*x2 + dd_recip_int[14]*x4 + dd_recip_int[15]*x6)
    z2 = x16*(dd_recip_int[17] + dd_recip_int[19]*x2 + dd_recip_int[21]*x4 + dd_recip_int[23]*x6)
    z3 = x24*(dd_recip_int[25] + dd_recip_int[27]*x2 + dd_recip_int[29]*x4 + dd_recip_int[31]*x6)

    (((z3+z2)+z1)+z+one(DD))*x
end

const dd_asinh_coeff = DD[
    DD(1.0),
    DD(-0.16666666666666666, -9.25185853854297e-18),
    DD(0.075, 2.7755575615628915e-18),
    DD(-0.044642857142857144, 9.912705577010326e-19),
    DD(0.030381944444444444, 3.854941057726238e-19),
    DD(-0.022372159090909092, 9.462128050782583e-19),
    DD(0.017352764423076924, -8.006416042969879e-19),
    DD(-0.01396484375, 6.938893903907229e-19),
    DD(0.011551800896139705, 8.163404592832033e-19),
    DD(-0.009761609529194078, -5.478074134663601e-19),
    DD(0.008390335809616815, 4.130293990420969e-19),
    DD(-0.0073125258735988454, 3.394024192128536e-19),
    DD(0.006447210311889649, -3.1225022567582527e-19),
    DD(-0.005740037670841924, 1.2849803525754126e-19),
    DD(0.005153309682319905, -3.888173308223878e-19),
    DD(-0.004660143486915096, 2.7979410902851727e-20),
    DD(0.004240907093679363, -1.5770213417970974e-19),
    DD(-0.003880964558837669, -1.858632295689436e-19),
    DD(0.0035692053938259347, -1.992587776459846e-19)
];

# |a| <= 0.045
function asinh_series(a::DD)
    x = a
    x2 = x*x
    x4 = x2*x2
    x6 = x2*x4
    x8 = x4*x4
    x16 = x8*x8
    x24 = x8*x16

    z = dd_asinh_coeff[2]*x2 + dd_asinh_coeff[3]*x4 + dd_asinh_coeff[4]*x6
    z1 = x8*(dd_asinh_coeff[5] + dd_asinh_coeff[6]*x2 + dd_asinh_coeff[7]*x4 + dd_asinh_coeff[8]*x6)
    z2 = x16*(dd_asinh_coeff[9] + dd_asinh_coeff[10]*x2 + dd_asinh_coeff[11]*x4 + dd_asinh_coeff[12]*x6)
    z3 = x24*(dd_asinh_coeff[13] + dd_asinh_coeff[14]*x2 + dd_asinh_coeff[15]*x4 + dd_asinh_coeff[16]*x6)

    (((z3+z2)+z1)+z+one(DD))*x
end


function atanh(x::DD)
    isneg, abs_x = signbit(x.hi), abs(x)
    if abs_x.hi > one(DD)
        throw(ErrorException("atanh arg $(x) outside -1..+1"))
        return dd_NaN
    elseif abs_x >= (one(DD)-1.0e-96)
        return isneg ? -dd_Inf : dd_Inf
    end
    if abs_x.hi <= 0.00275
        z = atanh_series(abs_x)
    else
        nm = one(TD) + abs_x
        dn = one(TD) - abs_x
        q = nm/dn
        q = log(DD(q))
        z = divby2(q)
    end
    isneg ? -z : z
end

function asinh(x::DD)
    isneg, abs_x = signbit(x),abs(x)
    if isneg
        -asinh(abs_x)
    elseif abs_x.hi <= 0.045
        b = asinh_series(abs_x)
    else
        y = one(TD)+sqrt(TD(x)*x+one(TD))
        y = TD(x) / y
        y = atanh(DD(y))
        mulby2(y)
    end
end

function acosh(x::DD)
   isneg, abs_x = signbit(x.hi), abs(x)
    if isneg || abs_x.hi < one(Float64)
        throw(ErrorException("acosh arg $(x) outside +1.."))
        return dd_NaN
    elseif abs_x <= (one(DD)+1.0e-96)
        return zero(DD)
    end

    if abs_x.hi <= 1.01 #1.00275
        a = sqrt(TD(abs_x)*abs_x - one(TD))
        b = asinh(DD(a))
        b
    else
        a = sqrt(TD(abs_x)*abs_x - one(TD))
        b = a + abs_x
        b = log(DD(b))
        b
    end
end

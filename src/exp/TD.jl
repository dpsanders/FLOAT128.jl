#=
#    Internal Use Only
=#

const n_inv_fact_TD = 32;

const inv_fact_TD = TD[
    #TD(1.0,0.0,0.0),
    #TD(0.5,0.0,0.0),
    TD(0.16666666666666666,9.25185853854297e-18,5.135813185032629e-34),
    TD(0.041666666666666664,2.3129646346357427e-18,1.2839532962581572e-34),
    TD(0.008333333333333333,1.1564823173178714e-19,1.6049416203226965e-36),
    TD(0.001388888888888889,-5.300543954373577e-20,-1.7386867553495878e-36),
    TD(0.0001984126984126984,1.7209558293420705e-22,1.4926912391394127e-40),
    TD(2.48015873015873e-5,2.1511947866775882e-23,1.865864048924266e-41),
    TD(2.7557319223985893e-6,-1.858393274046472e-22,8.491754604881993e-39),
    TD(2.755731922398589e-7,2.3767714622250297e-23,-3.263188903340883e-40),
    TD(2.505210838544172e-8,-1.448814070935912e-24,2.0426735146714455e-41),
    TD(2.08767569878681e-9,-1.20734505911326e-25,1.702227928892871e-42),
    TD(1.6059043836821613e-10,1.2585294588752098e-26,-5.31334602762985e-43),
    TD(1.1470745597729725e-11,2.0655512752830745e-28,6.889079232466646e-45),
    TD(7.647163731819816e-13,7.03872877733453e-30,-7.827539277162583e-48),
    TD(4.779477332387385e-14,4.399205485834081e-31,-4.892212048226615e-49),
    TD(2.8114572543455206e-15,1.6508842730861433e-31,-2.877771793074479e-50),
    TD(1.5619206968586225e-16,1.1910679660273754e-32,-4.577506059629983e-49),
    TD(8.22063524662433e-18,2.2141894119604265e-34,-1.508914023774199e-50),
    TD(4.110317623312165e-19,1.4412973378659527e-36,-5.285627548789812e-53),
    TD(1.9572941063391263e-20,-1.3643503830087908e-36,1.3392348251125064e-53),
    TD(8.896791392450574e-22,-7.911402614872376e-38,-3.1877976790570933e-54),
    TD(3.868170170630684e-23,-8.843177655482344e-40,3.8718157106173247e-56),
    TD(1.6117375710961184e-24,-3.6846573564509766e-41,1.613256546090552e-57),
    TD(6.446950284384474e-26,-1.9330404233703465e-42,-1.5213023807039144e-58),
    TD(2.4795962632247976e-27,-1.2953730964765229e-43,6.403390159849962e-60),
    TD(9.183689863795546e-29,1.4303150396787322e-45,-8.551226774650505e-62),
    TD(3.279889237069838e-30,1.5117542744029879e-46,8.058517719519716e-63),
    TD(1.1309962886447716e-31,1.0498015412959506e-47,-4.346150929397795e-64),
    TD(3.7699876288159054e-33,2.5870347832750324e-49,3.23789002742564e-66),
    TD(1.216125041553518e-34,5.586290567888806e-51,6.615948578082792e-68),
    TD(3.8003907548547434e-36,1.7457158024652518e-52,2.0674839306508725e-69)
];

function exp(a::TD)
  p2=16
  k = ldexp(1.0, p2)
  inv_k = 1.0 / k

  if (a.hi <= -709.0)
    return zero(Float64)
  elseif (a.hi >=  709.0)
    return td_Inf
  elseif (a.hi == zero(Float64))
    return one(TD)
  elseif (a == one(TD))
    return td_exp1
  end

  m = floor(a.hi / td_log2.hi + 0.5)
  #r = mul_pwr2(a - td_log2 * m, inv_k)
  r = a - td_log2*m
  r = ldexp(r, -p2)
  #TD s, p, t;
  # thresh = inv_k * TD::_eps;
  thresh = inv_k * eps(eps(1.0))

  p = sqr(r);
  s = r + divby2(p)
  i = 0
  p *= r
  t = p*inv_fact_TD[i+1]; i+=1
  s += t
  while (abs(t.hi) > thresh && i < 9)
    p *= r
    t = p*inv_fact_TD[i+1]; i+=1
    s += t
  end

  s = mulby2(s) + sqr(s);
  s = mulby2(s) + sqr(s);
  s = mulby2(s) + sqr(s);
  s = mulby2(s) + sqr(s);
  s = mulby2(s) + sqr(s);
  s = mulby2(s) + sqr(s);
  s = mulby2(s) + sqr(s);
  s = mulby2(s) + sqr(s);
  s = mulby2(s) + sqr(s);
  s = mulby2(s) + sqr(s);
  s = mulby2(s) + sqr(s);
  s = mulby2(s) + sqr(s);
  s = mulby2(s) + sqr(s);
  s = mulby2(s) + sqr(s);
  s = mulby2(s) + sqr(s);
  s = mulby2(s) + sqr(s);
  s += 1.0;

  m = trunc(m)

  ldexp(s, convert(Int,m));
end


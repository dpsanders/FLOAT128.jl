
function timeit(fn::Function, arg::Real, n::Int)
   mn = 9999999.9; md=0.0;
   a=fn(arg)
   for i in 1:n
      a = time_ns()
      b = fn(arg)
      c = time_ns()-a
      if (c < mn) mn = c end;
      md = md+c;
   end
   md1 = div(md,n)
   Int(mn),Int(md1),div(6*Int(mn)+Int(md1),7)
end


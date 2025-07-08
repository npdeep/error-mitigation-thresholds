function Oprod=ProductStateOp(Oi,M)


Oprod=1;

for m=1:M
    Oprod=kron(Oprod,sparse(Oi{m}));
end
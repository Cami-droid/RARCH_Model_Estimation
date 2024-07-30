function std_errors_calc(model,specification, outputs)

m_t=(mS_t',mD_t')'

mS_t=thetaS-vech(rt*rt')
mD_t=diff(log(ll(thetaS,thetaD),thetaD))

for RDCC model: 

m_rdcc_t=(mM_rdcc_t', mS_rdcc_t', mD_rdcc_t')'

mM_rdcc=zeros(1,d)

for j=1:d
mM_rdcc()=diff(log(ll(theta(j)),theta(j)))
end

mS_rdcc_t=thetaS_rdcc-vech8(et*et')

mD_rdcc_t=diff(log(ll(thetaM,thetaS,thetaD),thetaD))

J=var(1/sqrt(T))




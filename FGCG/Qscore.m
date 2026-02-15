function Q=Qscore(degreeset,A)

cov=sum(A,"all");
sp=1/sum(degreeset);
Q=cov*sp;

end
function zvec = pack_zvec(xval,uval,muval,lamval,gamval,slackval,T,N)
    zvec = [xval{1}];
    for t = 1:T
        zvec = [zvec; xval{t+1}];
        for i = 1:N
            zvec = [zvec; uval{t,i}; muval{t,i}; gamval{t,i}; slackval{t,i}; lamval{t,i}];
        end
    end
    for i = 1:N
        zvec = [zvec; muval{T+1,i}; gamval{T+1,i}; slackval{T+1,i}];
    end
end
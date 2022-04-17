function [ ] = DDE(NP,G,D)
    MAXRUN=20;
    global fes;
    fes=0;
    ok=0;
    % NP=20;
    % D=2;
    % G=100;
    F=0.5;
    CR=0.1;
    Xs=100;
    Xx=-100;
    error=10^-6;
    bestGen=zeros(1,MAXRUN);
    bestFes=zeros(1,MAXRUN);
    bestTime=zeros(1,MAXRUN);
    bestFun=zeros(1,MAXRUN);
    maxFes=10^5;
    if isempty(dir('DdeOutput'))
    else
        rmdir('DdeOutput', 's');
    end
    mkdir('DdeOutput');

    for run = 1:MAXRUN
        tic
        fes=0;
        x = zeros(D,NP);
        v = zeros(D,NP);
        u = zeros(D,NP);
        x = randi([Xx,Xs],D,NP);
    
        for m = 1:NP
            Ob(m) = func3(x(:,m));
        end
        trace(1) = max(Ob);
    
        fprintf('第%d次运行\t第%d代\t%g\n',run,0,trace(1));
        onerunfile = fopen(['.\DdeOutput\F1_run' num2str(run) '.txt'],'w'); %建立第run次运行的记录文档
        if fes==NP
            fprintf(onerunfile,'代数\t函数运行次数\t计算最大值\n');
        end
        fprintf(onerunfile,'%-d\t%10d\t\t%g\r\n',0,fes,trace(1));
        
        for gen = 1:G
            for m = 1:NP
                r1 = randi([1,NP],1,1);
                while r1 == m
                    r1 = randi([1,NP],1,1);
                end
                r2 = randi([1,NP],1,1);
                while r2 == m || r2 == r1
                    r2 = randi([1,NP],1,1);
                end
                r3 = randi([1,NP],1,1);
                while r3 == m || r3 == r1 || r3 == r2
                    r3 = randi([1,NP],1,1);
                end
                v(:,m) = floor(x(:,r1)+F*(x(:,r2)-x(:,r3)));
            end
            
        
            r = randi([1,D],1,1);
            for n = 1:D
                cr = rand(1);
                if (cr <= CR) | (n==r)
                    u(n,:) = v(n,:);
                else
                    u(n,:) = x(n,:);
                end
            end
        
            for n = 1:D
                for m = 1:NP
                    if u(n,m) < Xx
                        u(n,m) = Xx;
                    end
                    if u(n,m) > Xs
                        u(n,m) = Xs;
                    end
                end
            end
        
            for m = 1:NP
                Ob1(m) = func3(u(:,m));
            end
            for m = 1:NP
                if Ob1(m) > Ob(m)
                    x(:,m) = u(:,m);
                end
            end
        
            for m = 1:NP
                Ob(m) = func3(x(:,m));
            end
            trace(gen+1) = max(Ob);

            if fes >= maxFes
                fprintf(onerunfile,'%d\t%10d\t\t%g\r\n',gen,fes,trace(gen+1));
                break;
            end
            if (trace(gen)~=trace(gen+1) & abs(trace(gen)-trace(gen+1))<error)
                ok=ok+1;
                break;
            end
            if trace(gen+1) > trace(gen)
                bestGen(run) = gen;
                bestFes(run) = fes;
                bestFun(run) = trace(gen+1);
                bestTime(run) = toc;
            end
            fprintf(onerunfile,'%d\t%10d\t\t%g\r\n',gen,fes,trace(gen+1));
        end
        fclose(onerunfile);
        [SortOb,Index] = sort(Ob);
        X = x(:,Index);
        Xbest = X(:,end);
        Y = max(Ob);
    end

    F1total = fopen('.\DdeOutput\F1total.txt','w');
    for i = 1:MAXRUN
        fprintf(F1total,'第%d代:',i);
        fprintf(F1total,'(首次最优解%g)---',bestFun(i));
        fprintf(F1total,'(迭代数%d)---',bestGen(i));
        fprintf(F1total,'(函数评价次数%d)---',bestFes(i));
        fprintf(F1total,'(用时%gs)\r\n',bestTime(i));
    end
    fclose(F1total);

    bestfun=max(bestFun);
    worstfun=min(bestFun);
    meanfun=mean(bestFun);
    stdfun=std(bestFun,0);
    meanbestGen=mean(bestGen);
    meanbestFes=mean(bestFes);
    meanbestTime=mean(bestTime);
    okPercent=ok/MAXRUN*100;
    
    Total = fopen('.\DdeOutput\Total.txt','w');
    fprintf(Total,'测试函数名:DDE\r\n');
    fprintf(Total,'D:%d\r\n',D);
    fprintf(Total,'Pc:%.1f\r\n',CR);
    fprintf(Total,'Pm:%.1f\r\n',F);
    fprintf(Total,'NP:%d\r\n',NP);
    fprintf(Total,'MAXG:%d\r\n',G);
    fprintf(Total,'MAXRUN:%d\r\n',MAXRUN);
    fprintf(Total,'bestfun:%d\r\n',bestfun);
    fprintf(Total,'worstfun:%d\r\n',worstfun);
    fprintf(Total,'meanfun:%d\r\n',meanfun);
    fprintf(Total,'stdfun:%d\r\n',stdfun);
    fprintf(Total,'meanbestGen:%g\r\n',meanbestGen);
    fprintf(Total,'meanbestFes:%g\r\n',meanbestFes);
    fprintf(Total,'meanbestTime:%g\r\n',meanbestTime);
    fprintf(Total,'符合误差百分比(ok):%.2f%%\r\n',okPercent);
    fclose(Total);
end

function y = func3(x)
    y = -((x(1).^2+x(2)-1).^2+(x(1)+x(2).^2-7).^2)/200+10;
    global fes;
    fes = fes+1;
end
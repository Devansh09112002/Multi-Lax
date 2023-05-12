function [A] = combinator(N,K,s1,s2)

ng = nargin;

if ng == 2
    s1 = 'p';
    s2 = 'r';
elseif ng == 3 
    s2 = 'n';
elseif ng ~= 4
    error('Only 2, 3 or 4 inputs are allowed.  See help.')
end

if isempty(N) || K == 0
   A = [];  
   return
elseif numel(N)~=1 || N<=0 || ~isreal(N) || floor(N) ~= N 
    error('N should be one real, positive integer. See help.')
elseif numel(K)~=1 || K<0 || ~isreal(K) || floor(K) ~= K
    error('K should be one real non-negative integer. See help.')
end

STR = lower(s1(1)); % We are only interested in the first letter.

if ~strcmpi(s2(1),'r')
    STR = [STR,'n'];
else
   STR = [STR,'r']; 
end

try
    switch STR
        case 'pr'
            A = perms_rep(N,K);     % strings
        case 'pn'
            A = perms_no_rep(N,K);  % permutations
        case 'cr'
            A = combs_rep(N,K);     % multichoose
        case 'cn'
            A = combs_no_rep(N,K);  % choose
        otherwise
            error('Unknown option passed.  See help')
    end
catch
    rethrow(lasterror) % Throw error from here, not subfunction.
    % The only error thrown should be K>N for non-replacement calls.
end

function PR = perms_rep(N,K)
        
%      tic,Tc = combinator(single(9),7,'p','r');toc  
%      %Elapsed time is 0.199397 seconds.  Allow Ctrl+T+C+R on block
%      tic,Tj = combn(single(1:9),7);toc  
%      %Elapsed time is 0.934780 seconds.
%      isequal(Tc,Tj)  % Yes

if N==1
   PR = ones(1,K,class(N)); 
   return
elseif K==1
    PR = (1:N).';
    return
end

CN = class(N);
M = double(N);  % Single will give us trouble on indexing.
L = M^K;  % This is the number of rows the outputs will have.
PR = zeros(L,K,CN);  % Preallocation.
D = ones(1,N-1,CN);  % Use this for cumsumming later.
LD = M-1;  % See comment on N. 
VL = [-(N-1) D].';  % These values will be put into PR.
    
% Now start building the matrix.
TMP = VL(:,ones(L/M,1,CN));  % Instead of repmatting.
PR(:,K) = TMP(:);  % We don't need to do two these in loop.
PR(1:M^(K-1):L,1) = VL;  % The first column is the simplest.
% Here we have to build the cols of PR the rest of the way.
for ii = K-1:-1:2
    ROWS = 1:M^(ii-1):L;  % Indices into the rows for this col.
    TMP = VL(:,ones(length(ROWS)/(LD+1),1,CN));  % Match dimension.
    PR(ROWS,K-ii+1) = TMP(:);  % Build it up, insert values.
end

PR(1,:) = 1;  % For proper cumsumming.
PR = cumsum2(PR);  % This is the time hog.

function PN = perms_no_rep(N,K)
        
% Subfunction: permutations without replacement.
% Uses the algorithm in combs_no_rep as a basis, then permutes each row.
% pn = @(N,K) prod(1:N)/(prod(1:(N-K)));  Number of rows.

if N==K
    PN = perms_loop(N);  % Call helper function.
%     [id,id] = sort(PN(:,1));  %#ok  Not nec., uncomment for nice order.
%     PN = PN(id,:);  % Return values.
    return
elseif K==1
    PN = (1:N).';  % Easy case.
    return
end

if K>N  % Since there is no replacement, this cannot happen.
    error(['When no repetitions are allowed, '...
           'K must be less than or equal to N'])
end

M = double(N);  % Single will give us trouble on indexing.
WV = 1:K;  % Working vector.
lim = K;   % Sets the limit for working index.
inc = 1;   % Controls which element of WV is being worked on.
BC = prod(M-K+1:M);  % Pre-allocation of return arg.
BC1 = BC / ( prod(1:K)); % Number of comb blocks.
PN = zeros(round(BC),K,class(N));
L = prod(1:K) ;  % To get the size of the blocks.
cnt = 1+L;
P = perms_loop(K);  % Only need to use this once.
PN(1:(1+L-1),:) = WV(P);  % The first row.

for ii = 2:(BC1 - 1);
    if logical((inc+lim)-N)  % The logical is nec. for class single(?)
        stp = inc;  % This is where the for loop below stops.
        flg = 0;  % Used for resetting inc.
    else
        stp = 1;
        flg = 1;
    end
    
    for jj = 1:stp
        WV(K  + jj - inc) = lim + jj;  % Faster than a vector assignment!
    end
                                                                             
    PN(cnt:(cnt+L-1),:) = WV(P);  % Assign block.
    cnt = cnt + L;  % Increment base index.    
    inc = inc*flg + 1;  % Increment the counter.
    lim = WV(K - inc + 1 );  % lim for next run.
end

V = (N-K+1):N;  % Final vector.
PN(cnt:(cnt+L-1),:) = V(P);  % Fill final block.

function P = perms_loop(N)
%
%      tic,Tp = perms(1:9);toc
%      %Elapsed time is 0.222111 seconds.  Allow Ctrl+T+C+R on block
%      tic,Tc = combinator(9,9,'p');toc  
%      %Elapsed time is 0.143219 seconds.
%      isequal(Tc,Tp)  % Yes

M = double(N); % Single will give us trouble on indexing.
P = 1;  % Initializer.
G = cumprod(1:(M-1));  % Holds the sizes of P.
CN = class(N);

for n = 2:M
    q = P;
    m = G(n-1);
    P = zeros(n*m,n,CN);
    P(1:m, 1) = n;
    P(1:m, 2:n) = q;
    a = m + 1;
    
    for ii = n-1:-1:1,
        t = q;
        t(t == ii) = n;
        b = a + m - 1;
        P(a:b, 1) = ii;
        P(a:b, 2:n) = t;
        a = b + 1;
    end 
end

function CR = combs_rep(N,K)
% Subfunction multichoose:  combinations with replacement.
% cr = @(N,K) prod((N):(N+K-1))/(prod(1:K)); Number of rows.

M = double(N);  % Single will give us trouble on indexing.
WV = ones(1,K,class(N));  % This is the working vector.
mch = prod((M:(M+K-1)) ./ (1:K));  % Pre-allocation.
CR = ones(round(mch),K,class(N));

for ii = 2:mch
    if WV(K) == N
        cnt = K-1;  % Work backwards in WV.
        
        while WV(cnt) == N
            cnt = cnt-1;  % Work backwards in WV.
        end

        WV(cnt:K) = WV(cnt) + 1;  % Fill forward.
    else
        WV(K) = WV(K)+1;   % Keep working in this group.
    end

    CR(ii,:) = WV;
end

function CN = combs_no_rep(N,K)
    
% cn = @(N,K) prod(N-K+1:N)/(prod(1:K));  Number of rows.
% Same output as the MATLAB function nchoosek(1:N,K), but often faster for
% larger N.
%
%      tic,Tn = nchoosek(1:17,8);toc
%      %Elapsed time is 0.430216 seconds.  Allow Ctrl+T+C+R on block
%      tic,Tc = combinator(17,8,'c');toc  
%      %Elapsed time is 0.024438 seconds.
%      isequal(Tc,Tn)  % Yes

if K>N
    error(['When no repetitions are allowed, '...
           'K must be less than or equal to N'])
end

M = double(N);  % Single will give us trouble on indexing.

if K == 1
   CN =(1:N).';  % These are simple cases.
   return
elseif K == N
    CN = (1:N);
    return
elseif K==2 && N>2  % This is an easy case to do quickly.
    BC = (M-1)*M / 2;
    id1 = cumsum2((M-1):-1:2)+1;
    CN = zeros(BC,2,class(N));
    CN(:,2) = 1;
    CN(1,:) = [1 2];
    CN(id1,1) = 1;
    CN(id1,2) = -((N-3):-1:0);
    CN = cumsum2(CN);
    return
end

WV = 1:K;  % Working vector.
lim = K;   % Sets the limit for working index.
inc = 1;   % Controls which element of WV is being worked on.
BC = prod(M-K+1:M) / (prod(1:K));  % Pre-allocation.
CN = zeros(round(BC),K,class(N));
CN(1,:) = WV;  % The first row.

for ii = 2:(BC - 1);   
    if logical((inc+lim)-N) % The logical is nec. for class single(?)
        stp = inc;  % This is where the for loop below stops.
        flg = 0;  % Used for resetting inc.
    else
        stp = 1;
        flg = 1;
    end
    
    for jj = 1:stp
        WV(K  + jj - inc) = lim + jj;  % Faster than a vector assignment.
    end
    
    CN(ii,:) = WV;  % Make assignment.
    inc = inc*flg + 1;  % Increment the counter.
    lim = WV(K - inc + 1 );  % lim for next run. 
end
  
CN(ii+1,:) = (N-K+1):N;

function A = cumsum2(A)

if isfloat(A)
    A = cumsum(A);  % For single and double, use built-in.
    return
else
    try
        A = cumsumall(A);  % User has the MEX-File ready?
    catch
        warning('Cumsumming by loop.  MEX cumsumall.cpp for speed.') %#ok
        for ii = 2:size(A,1)
            A(ii,:) = A(ii,:) + A(ii-1,:); % User likes it slow.
        end
    end
end

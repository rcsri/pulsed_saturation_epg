function [Mz_out, Mxy_out] = bmsim_sequence_epg_ver2(p0, sequence, lstype)
% Simulate Bloch McConnell exchange pulse sequence with extended phase graphs

fsat            = sequence.fsat;
t_sat           = sequence.t_sat;

b1_sat          = sequence.b1_sat; % shaped pulse, scaled to 1
if isfield(sequence,'b1_power')
    b1_power          = sequence.b1_power; % peak power uT
else
    % this implies that 'b1_sat' is scaled to the right power,
    % so doesn't need an extra 'b1_power' scaling factor
    % however, if b1_sat is scaled to 1 (e.g. from waveform file),
    % then use b1_power to scale back to appropriate uT.
    b1_power = 1;
    sequence.b1_power = 1;
end

gamma_          = sequence.gamma_;
num_sat         = sequence.num_sat;
tspoil_sat      = sequence.tspoil_sat;

if strcmp(sequence.readout_type, 'TFE')
    read_seq = sequence.tfe;
elseif strcmp(sequence.readout_type, 'EPI')
    read_seq = sequence.epi;
end

tex                 = read_seq.tex;
num_ex              = read_seq.num_ex;
num_ex_dummy        = read_seq.num_ex_dummy;
flipex_vect_dummy   = read_seq.flipex_vect_dummy;
flipex_vect_ex      = read_seq.flipex_vect_ex;
tspoil_tfe          = read_seq.tspoil;
tdel_pre            = sequence.tdel_pre;
tdel_end            = sequence.tdel_end;

% voxeldependent
b1_scale        = sequence.b1_scale;    % additional scale factor for b1 map, for saturation
b1_scale_ex     = sequence.b1_scale;  %_ex; % additional scale factor for b1 map, for excitation
b0              = sequence.b0;

% validate input
if ~is_valid_pool_params(p0)
    disp('invalid pool parameters');
end

% b0 correction
p0(:,1) = p0(:,1) + b0;

nf = numel(fsat);

% Extract equilibrium and inhomogeneous term
[Meq, Ceq] = bmsim_Meq( p0 );
I = eye(numel(Meq));

Mz_out = zeros(1,nf);
Mxy_out = zeros(1,nf);

noadd = 1;
Nstates = 25;

Ns = 5; % max Ns = Nstates, Ns < Nstates to save time.

num_compounds   = size(p0,1);
num_MT          = p0(end,6);
num_components  = 3 * num_compounds - (num_MT * 2);  % number of rows in A

N = num_compounds - num_MT; % water or cest pools

for ixc = 1:num_compounds
    is_compound_MT = (p0(ixc,6) == 1);
    U0       =  [ 1 +1i 0;  1  -1i   0;  0  0 1];
    U = kron(eye(N), U0);
    if is_compound_MT
        U = blkdiag( U, 1 );
    end
    inv_U = inv(U);
end
FZ      = zeros(num_components,Nstates);	    % State matrix
FZ(:,1) = Meq;

C = zeros(num_components,Nstates);
C(:,1) = Ceq;


%----------------------------------------------------------------------%
% find repeated freqs to avoid repeatedly calculating the same matrices
%----------------------------------------------------------------------%
% Find the unique values of A and the index vectors ia and ic, such that C = A(ia) and A = C(ic).
% [C, ia, ic] = unique(A)

f_unique = unique(fsat);
if (length(f_unique) < length(fsat))
    has_repeated_f = 1;
else
    has_repeated_f = 0;
end
if has_repeated_f
    DICT_expm_UAU_sat       = cell(length(f_unique),1);
    DICT_inv_UAU_sat        = cell(length(f_unique),1);
    DICT_expm_UAU_sat_relax = cell(length(f_unique),1);
    DICT_inv_UAU_sat_relax  = cell(length(f_unique),1);
    DICT_expm_UAU_read_relax = cell(length(f_unique),1);
    DICT_inv_UAU_read_relax = cell(length(f_unique),1);
    DICT_expm_UAU_pre       = cell(length(f_unique),1);
    DICT_inv_UAU_pre       = cell(length(f_unique),1);
    DICT_expm_UAU_end       = cell(length(f_unique),1);
    DICT_inv_UAU_end       = cell(length(f_unique),1);
    DICT_E = cell(length(f_unique),1);
    DICT_B_sat = cell(length(f_unique),1);
    DICT_B_rel = cell(length(f_unique),1);
end

%----------------------------------------------%
% loop through offsets
%-----------------------------------------------%
for ix_f = 1:nf
    % reset magnetization
    
    f = fsat(ix_f);
    
    %===================================================
    % PREP
    %===================================================
    
    % sat propagators
    w1_sat   = gamma_ * b1_scale * b1_power * b1_sat;  % note: b1_sat can be a vector if shaped pulse
    
    ind_dict = find(f_unique == f);
    val = DICT_expm_UAU_sat{ind_dict};
    
    if isempty(val) % if does not exist, store in dictionary
        [expm_UAU_sat, inv_UAU_sat] = bmsim_compute_prop_shaped_pulse_epg(p0, f, w1_sat, t_sat, U, inv_U, lstype);
        
        if has_repeated_f
            DICT_expm_UAU_sat{ind_dict} = expm_UAU_sat;
            DICT_inv_UAU_sat{ind_dict} = inv_UAU_sat;
        end
    else % if already exists, retrieve from dictionary
        expm_UAU_sat = DICT_expm_UAU_sat{ind_dict};
        inv_UAU_sat = DICT_inv_UAU_sat{ind_dict};
    end
    
    % relaxation propagators
    if isempty(val)
        A_relax                         = bmsim_mtx( p0, f, 0 , lstype);
        UAU_relax                       = U * A_relax * inv_U;
        [expm_UAU_sat_relax, inv_UAU_sat_relax]     = bmsim_compute_prop( tspoil_sat,   UAU_relax);        % relax in between sat
        [expm_UAU_read_relax, inv_UAU_read_relax]   = bmsim_compute_prop( tspoil_tfe,   UAU_relax);        % relax & spoil for TFE
        [expm_UAU_pre, inv_UAU_pre]                 = bmsim_compute_prop( tdel_pre,     UAU_relax);        % relax pre-readout
        [expm_UAU_end, inv_UAU_end]                 = bmsim_compute_prop( tdel_end,     UAU_relax);        % relax post-readout
        
        if has_repeated_f
            DICT_expm_UAU_sat_relax{ind_dict}     = expm_UAU_sat_relax;
            DICT_inv_UAU_sat_relax{ind_dict}      = inv_UAU_sat_relax;
            DICT_expm_UAU_read_relax{ind_dict}    = expm_UAU_read_relax;
            DICT_inv_UAU_read_relax{ind_dict}     = inv_UAU_read_relax;
            DICT_expm_UAU_pre{ind_dict}           = expm_UAU_pre;
            DICT_inv_UAU_pre{ind_dict}            = inv_UAU_pre;
            DICT_expm_UAU_end{ind_dict}           = expm_UAU_end;
            DICT_inv_UAU_end{ind_dict}            = inv_UAU_end;
        end
    else
        expm_UAU_sat_relax          = DICT_expm_UAU_sat_relax{ind_dict};
        inv_UAU_sat_relax           = DICT_inv_UAU_sat_relax{ind_dict};
        expm_UAU_read_relax         = DICT_expm_UAU_read_relax{ind_dict};
        inv_UAU_read_relax          = DICT_inv_UAU_read_relax{ind_dict};
        expm_UAU_pre                = DICT_expm_UAU_pre{ind_dict};
        inv_UAU_pre                 = DICT_inv_UAU_pre{ind_dict};
        expm_UAU_end                = DICT_expm_UAU_end{ind_dict};
        inv_UAU_end                 = DICT_inv_UAU_end{ind_dict};
    end
    
    
    %===================================================
    % PROPAGATE
    %===================================================
    
    FZ(:,2:end) = 0; % keep only F0 and Z0 from previous TR
    
    %---------------------------------------------------------------------%
    % SATURATION (EPG)
    %---------------------------------------------------------------------%
    num_intervals = size(expm_UAU_sat, 3);
    if isempty(val)
        if num_intervals > 1
            E       = zeros(size(expm_UAU_sat,1),size(expm_UAU_sat,2), num_intervals);
            B_sat   = zeros(size(expm_UAU_sat,1), num_intervals);
            for ind = 1:num_intervals
                E(:,:,ind)      = expm_UAU_sat_relax * expm_UAU_sat(:,:,ind);
                B_sat(:,ind)    = inv_UAU_sat(:,:,ind) * ( expm_UAU_sat(:,:,ind) - I) * C(:,1);
            end
            B_rel = inv_UAU_sat_relax * ( expm_UAU_sat_relax - I) * C(:,1);
        else
            E = expm_UAU_sat_relax * expm_UAU_sat;
            B_sat = inv_UAU_sat * ( expm_UAU_sat - I) * C(:,1);
            B_rel = inv_UAU_sat_relax * ( expm_UAU_sat_relax - I) * C(:,1);
        end
        
        DICT_E{ind_dict}        = E;
        DICT_B_sat{ind_dict}    = B_sat;
        DICT_B_rel{ind_dict}    = B_rel;
        
    else
        E       = DICT_E{ind_dict};
        B_sat   = DICT_B_sat{ind_dict};
        B_rel   = DICT_B_rel{ind_dict};
    end
    
    
    for ix_sat = 1:num_sat
        for ind = 1:num_intervals
            FZ(:,1) = bmsim_propagate_ver2(FZ(:,1), expm_UAU_sat(:,:,ind), B_sat(:,ind) );
            FZ(:,2:Ns) = E(:,:,ind) * FZ(:,2:Ns);
        end
        FZ(:,1) = bmsim_propagate_ver2(FZ(:,1), expm_UAU_sat_relax, B_rel );
        FZ(1:3,:) = epg_grad(FZ(1:3,:),noadd);
    end
    
    FZ(:,2:end) = 0; % keep only F0 and Z0
    
    %-----------------------%
    % delay, pre-readout
    %-----------------------%
    FZ(:,1) = bmsim_propagate(FZ(:,1), expm_UAU_pre, inv_UAU_pre, C(:,1));
    FZ(:,2:Ns) = expm_UAU_pre * FZ(:,2:Ns); % s>1 where C = 0, combine sat and relax into E
    
    %-----------------------%
    % dummy ex
    %-----------------------%
    for ix_ex = 1:num_ex_dummy
        flipex                          = flipex_vect_dummy(ix_ex);
        A_ex                            = bmsim_mtx( p0, 0, b1_scale_ex * deg2rad(flipex)/tex, lstype);
        UAU_ex                          = U * A_ex * inv_U;
        [expm_UAU_ex, inv_UAU_ex]       = bmsim_compute_prop( tex, UAU_ex);
        
        FZ(:,1) = bmsim_propagate(FZ(:,1), expm_UAU_ex, inv_UAU_ex, C(:,1));
        FZ(:,1) = bmsim_propagate(FZ(:,1), expm_UAU_read_relax, inv_UAU_read_relax, C(:,1));
        E = expm_UAU_read_relax * expm_UAU_ex;
        FZ(:,2:Ns) = E * FZ(:,2:Ns);
        FZ(1:3,:) = epg_grad(FZ(1:3,:),noadd); % spoil
    end
    
    %-----------------------%
    % readout ex
    %-----------------------%
    for ix_ex = 1:num_ex
        
        flipex = flipex_vect_ex(ix_ex);
        A_ex                            = bmsim_mtx( p0, 0, b1_scale_ex * deg2rad(flipex)/tex, lstype);
        UAU_ex                          = U * A_ex * inv_U;
        [expm_UAU_ex, inv_UAU_ex]       = bmsim_compute_prop( tex, UAU_ex);
        FZ(:,1) = bmsim_propagate(FZ(:,1), expm_UAU_ex, inv_UAU_ex, C(:,1));
        FZ(:,2:Ns) = expm_UAU_ex * FZ(:,2:Ns);
        
        %----------------------%
        % record magnetization
        %----------------------%
        if ix_ex == 1  % centric TFE
            M = inv_U * FZ(:,1);
            Mz_out(ix_f)    = real(M(3));
            Mxy_out(ix_f)   = abs(M(2));
        end
        
        FZ(:,1) = bmsim_propagate(FZ(:,1), expm_UAU_read_relax, inv_UAU_read_relax, C(:,1));
        FZ(:,2:Ns) = expm_UAU_read_relax * FZ(:,2:Ns);
        FZ(1:3,:) = epg_grad(FZ(1:3,:),noadd);
        
    end  % ix_ex
    
    % delay after readout
    FZ(:,1) = bmsim_propagate(FZ(:,1), expm_UAU_end, inv_UAU_end, C(:,1));
    FZ(:,2:Ns) = expm_UAU_end * FZ(:,2:Ns);
        
end % ix_f


function [FinalIndices, VarIdx, ExpDesign, Bstr, Cost] = uq_sobol_greater_order(Order, Sample1, Sample2, M_Sample1, Options, Bstr, CurrentAnalysis)
% [INDICES,VARIDX,ED,BSTR,COST] = UQ_SOBOL_GREATER_ORDER(ORDER, SAMPLE1, 
%                   SAMPLE2, M_SAMPLE1, OPTIONS, BSTR, ANALYSISOBJ):
%     calculate sample-based Sobol' indices of order ORDER for the analysis
%     object ANALYSISOBJ based on the samples SAMPLE1,SAMPLE2 and mixed
%     sample M_SAMPLE1. Additional options are specified in OPTIONS,
%     bootstrap options in BSTR.
%
% See also: UQ_SOBOL_INDEX,UQ_SOBOL_INDICES

% Size of the samples:
N = size(Sample1, 1);

% Number of factors:
FactorIndex = Options.FactorIndex;
M = length(FactorIndex);
NFactors = sum(FactorIndex);

VarIdxFull = nchoosek(1:M, Order);
NofIndicesFull = size(VarIdxFull, 1);

% Check the indices we need:
VarIdx = nchoosek(find(FactorIndex), Order);
NofIndices = size(VarIdx, 1);

% Number of samples that we have to get out:
NofSamples = size(VarIdx, 1);

% Bootstrap replications:
B = Options.Bootstrap.Replications;

% Initialize bootstrap variables:
if B > 0
    % Auxiliar matrix for bootstrap:
    AuxMatrix = [];
    for i = 1:NofSamples
        AuxMatrix = [AuxMatrix (i-1)*N*ones(1,N)];
    end
    
    % Remove output from the indices functions during Bootstrap:
    PseudoOptions = Options;
    PseudoOptions.Display = 0;
    
    % Allocate memory for the bootstrap estimates:
    B_Indices = zeros(B, NofSamples);
end


% First, create the sample needed, they are based on Sample2:
Shuffled_Sample = repmat(Sample2, NofSamples, 1);

% Now, start freezing the values needed:
Count = 0;
RealPos = [];
for ii = 1:NofIndicesFull
    % Select the parts of the sample that should be fixed
    VectorIdx = VarIdxFull(ii,:);
    SkipIndex = false;
    for jj = 1:length(FactorIndex)
        if ~FactorIndex(jj) && any(jj == VectorIdx)
            SkipIndex = true;
        end
    end
    
    if SkipIndex
        continue
        
    else
        % Bookkeeping to store the evaluations
        RealPos = [RealPos, (ii - 1)*N + 1:ii*N];
        Count = Count + 1;
    end
       
    % Replace them in the sample:
    for jj = 1:length(VectorIdx)
        Shuffled_Sample((Count - 1)*N + 1:Count*N, VectorIdx(jj)) = Sample1(:, VectorIdx(jj));
    end
end

M_Shuffled_Sample_Small = uq_evalModel(Options.Model,Shuffled_Sample);
Cost = size(M_Shuffled_Sample_Small, 1);

% Number of outputs:
NOuts = size(M_Shuffled_Sample_Small, 2);


% Retrieve the evaluations:
if Options.SaveEvaluations
    ExpDesign.X = Shuffled_Sample;
    ExpDesign.Y = M_Shuffled_Sample_Small;
else
    clear Shuffled_Sample
    ExpDesign = [];
end

M_Shuffled_Sample = zeros(NofIndicesFull*N, NOuts);
M_Shuffled_Sample(RealPos, :) = M_Shuffled_Sample_Small;

% Preallocate the outputs
TotalGroupIndices = zeros(NOuts, NofIndices);

for oo = 1:NOuts % Loop on the outputs:
    
    % Number of orders:
    NofOrders = length(CurrentAnalysis.Internal.Results.AllOrders);
    
    % Get the total indices for the group:
    tgindex = ...
        uq_sobol_index(Options, M_Sample1(:, oo), M_Shuffled_Sample(:, oo), ...
        CurrentAnalysis.Internal.f0(oo), CurrentAnalysis.Internal.D(oo), ...
        FactorIndex, VarIdxFull);
    
    % Due to "unused" variables NaNs will be returned. They need to be
    % removed and only the TotalGroup indices that have to do with the
    % used variables are saved:
    TotalGroupIndices(oo, :) = tgindex(~isnan(tgindex));
   
    % Check that all of them are there:
    if NofOrders ~= Order - 1
        error('Cannot calculate indices of order %d, previous indices are missing', ...
            Order);
    end
    
    

    % High order indices need lower order indices to be substracted. 
    SubstractIndex = cell(NofIndices, Order - 1);
    for jj = 1:NofIndices
        CurrentVars = VarIdx(jj,:);
        for ii = 1:NofOrders
            % ii refers to the previous indices variable
            PreviousVarIdx = CurrentAnalysis.Internal.Results.VarIdx{ii};
            PreviousIndex = CurrentAnalysis.Internal.Results.AllOrders{ii}(:,oo);
            
            % Loop inside the indices of this order:
            IdxToSubstract = [];
            for hh = 1 : size(PreviousVarIdx, 1)
                % Now we have to compare the current VarIdx with PreviousVarIdx
                % if all the PreviousVarIdx elements are in VarIdx, then this
                % index must be substracted
                
                % Vector containing the variables of this index
                PreviousVars = PreviousVarIdx(hh,:);
                
                % To avoid a loop for comparing them, let's repeat the vectors
                % along each other's dimension and check  that the resulting
                % matrix have the same number of ones as PreviousVars
                RepeatPrevious = repmat(PreviousVars', 1, length(CurrentVars));
                RepeatCurrent = repmat(CurrentVars, length(PreviousVars), 1);
                if sum(sum(RepeatPrevious == RepeatCurrent)) == size(PreviousVars, 2)
                    IdxToSubstract = [IdxToSubstract, hh];
                    TotalGroupIndices(oo, jj) = TotalGroupIndices(oo, jj) - PreviousIndex(hh);
                end
            end
            SubstractIndex{jj, ii} = IdxToSubstract;
        end
    end
    %% Bootstrap!
    if B > 0
        for bsample = 1:B
            % Reuse previous resampling:
            CurrentReSample = CurrentAnalysis.Internal.ReSampleIdx(bsample, :);
            B_Sample1 = M_Sample1(CurrentReSample, oo);
            
            % Map the results vectors to the new subsample:
            IdxReSample = repmat(CurrentReSample,1,NofSamples) + AuxMatrix;
            
            B_Shuffled_Sample = M_Shuffled_Sample(IdxReSample, oo);

            % Compute independent terms:
            Bf0 = sum(B_Sample1)/N;
            BD = (sum(B_Sample1.^2))/N - Bf0^2;

            % Bootstrap indices
            B_Indices(bsample, :) = ...
                uq_sobol_index(PseudoOptions, B_Sample1,...
                B_Shuffled_Sample, Bf0, BD, FactorIndex, VarIdx);
            
            % Subtract the previous indices:
            for jj = 1:NofIndices
                for ii = 1:Order - 1
                    ToSubstract = SubstractIndex{jj, ii};
                    if ~isempty(ToSubstract)
                        for hh = 1:length(ToSubstract)
                            B_Indices(bsample, jj) = ....
                                B_Indices(bsample, jj) - CurrentAnalysis.Internal.Bstr.AllOrders{ii}.Mean(ToSubstract(hh));
                        end
                    end
                end
            end
        end
        
        % Check for the confidence intervals and confidence level:
        [BGreater.CI, BGreater.Conf, BGreater.Mean] = ...
            uq_Bootstrap_CI(B_Indices, Options.Bootstrap.Alpha);
        
        % Assemble the bootstrap results
        Bstr.AllOrders{Order}.CI(:,oo,:) = BGreater.CI';
        Bstr.AllOrders{Order}.Mean(:,oo) = BGreater.Mean;
        Bstr.AllOrders{Order}.Conf(:,oo) = BGreater.Conf;
    else
        Bstr = [];
    end
    
end

%% Return value
FinalIndices = TotalGroupIndices;



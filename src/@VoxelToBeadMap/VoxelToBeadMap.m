classdef VoxelToBeadMap < handle
    properties
        Map (:,2) int32
    end

    methods (Access = public)
        function this = VoxelToBeadMap(beadIdxList)
            this.Map = [cat(1, beadIdxList{:}), ...
                        repelem(1 : numel(beadIdxList), cellfun(@(x) numel(x), beadIdxList))'];
            this.Map = sortrows(this.Map);
        end

        function [bIdx] = getBeadIdx(this, vIdx)
            bIdx = binarySearch(this.Map(:, 1), vIdx);
            if bIdx > 0
                bIdx = this.Map(bIdx, 2);
            end
        end
    end
end
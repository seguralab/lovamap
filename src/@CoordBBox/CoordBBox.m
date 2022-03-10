classdef CoordBBox
    properties
        Min (1,3) int32 = intmax('int32')
        Max (1,3) int32 = intmin('int32')
    end

    properties (Dependent, SetAccess = private)
        Dim    (1,3) uint32
        Volume (1,1) uint64
        IJKs   (:,3) int32
    end

    methods (Access = public)
        function this = CoordBBox(aMin, aMax)
            if nargin == 0
            elseif nargin == 1
                this.Min = aMin;
                this.Max = aMin;
            elseif nargin == 2
                this.Min = aMin;
                this.Max = aMax;
            else
                error('CoordBBox constructor requires 0, 1, or 2 arguments.');
            end
        end

        % Check if this bounding box is empty
        function b = isEmpty(this)
            b = any(this.Max < this.Min);
        end

        % Pad this bounding box with the specified padding
        function this = expand(this, aPad)
            this.Min = this.Min - int32(aPad);
            this.Max = this.Max + int32(aPad);
        end
        
        % Intersect this bounding box with the given bounding box.
        function this = intersect(this, aBBox)
            this.Min = max(this.Min, aBBox.Min);
            this.Max = min(this.Max, aBBox.Max);
        end
        
        % Translate this bounding box by the specified delta.
        function this = translate(this, aDelta)
            this.Min = this.Min + aDelta;
            this.Max = this.Max + aDelta;
        end
    end
    
    methods
        function vol = get.Volume(this)
            if this.isEmpty()
                vol = uint64(0);
            else
                vol = prod(uint64(this.Max - this.Min) + 1, 'native');
            end
        end
        
        function dim = get.Dim(this)
            if this.isEmpty()
                dim = zeros(1, 3, 'uint32');
            else
                dim = uint32(this.Max - this.Min) + 1;
            end
        end
        
        function ijks = get.IJKs(this)
            if this.isEmpty()
                ijks = [];
            else
                [ii, jj, kk] = ndgrid(this.Min(1) : this.Max(1), ...
                                      this.Min(2) : this.Max(2), ...
                                      this.Min(3) : this.Max(3));
                ijks = [ii(:) jj(:) kk(:)];
            end
        end
    end
end

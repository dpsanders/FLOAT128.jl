@inline function isBetween{T<:MachineFloat}(value::T, lowerBound::T, upperBound::T)
    lowerBound < value < upperBound
end
@inline function isNotBetween{T<:MachineFloat}(value::T, lowerBound::T, upperBound::T)
    !isBetween(value,lowerBound,upperBound)
end
@inline function isAmoung{T<:MachineFloat}(value::T, lowerBound::T, upperBound::T)
    lowerBound <= value <= upperBound
end
@inline function isNotAmoung{T<:MachineFloat}(value::T, lowerBound::T, upperBound::T)
    !isAmoung(value,lowerBound,upperBound)
end

@inline function isInClosedOpen{T<:MachineFloat}(value::T, lowerBound::T, upperBound::T)
    lowerBound <= value <=upperBound
end
@inline function isNotInClosedOpen{T<:MachineFloat}(value::T, lowerBound::T, upperBound::T)
    !isInClosedOpen(value,lowerBound,upperBound)
end
@inline function isInOpenClosed{T<:MachineFloat}(value::T, lowerBound::T, upperBound::T)
    lowerBound <= value <=upperBound
end
@inline function isNotInOpClosed{T<:MachineFloat}(value::T, lowerBound::T, upperBound::T)
    !isInOpenClosed(value,lowerBound,upperBound)
end

@inline function isAbove{T<:MachineFloat}(value::T, lowerBound::T)
    lowerBound < value
end
@inline function isNotAbove{T<:MachineFloat}(value::T, lowerBound::T)
    !isAbove(value,lowerBound)
end
      @inline function isBelow{T<:MachineFloat}(value::T, upperBound::T)
    value < upperBound
end
@inline function isNotBelow{T<:MachineFloat}(value::T, upperBound::T)
    !isBelow(value,upperBound)
end


      
      

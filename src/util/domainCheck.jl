immutable BotTop{T <: Real}
    bot::T
    top::T
    
    function BotTop(bot::T, top::T)
        bot, top = min(bot,top), max(bot,top)
        new(bot, top)
    end
end    

@inline function isBetween{T<:Real}(value::T, bounds::BotTop{T})
    bounds.bot < value < bounds.top
end
@inline function isNotBetween{T<:Real}(value::T, bounds::BotTop{T})
    !isBetween(value,bounds)
end
@inline function isAmoung{T<:Real}(value::T, bounds::BotTop{T})
    bounds.bot <= value <= bounds.top
end
@inline function isNotAmoung{T<:Real}(value::T, bounds::BotTop{T})
    !isAmoung(value,bounds)
end

@inline function isInClosedOpen{T<:Real}(value::T, bounds::BotTop{T})
    bounds.bot <= value < bounds.top
end
@inline function isNotInClosedOpen{T<:Real}(value::T, bounds::BotTop{T})
    !isInClosedOpen(value,lowerBound,upperBound)
end
@inline function isInOpenClosed{T<:Real}(value::T, bounds::BotTop{T})
    bounds.bot < value <= bounds.top
end
@inline function isNotInOpenClosed{T<:Real}(value::T, bounds::BotTop{T})
    !isInOpenClosed(value,lowerBound,upperBound)
end

@inline function isAbove{T<:Real}(value::T, bounds::BotTop{T})
    bounds.bot < value
end
@inline function isNotAbove{T<:Real}(value::T, bounds::BotTop{T})
    !isAbove(value,lowerBound)
end
@inline function isBelow{T<:Real}(value::T, bounds::BotTop{T})
    value < bounds.top
end
@inline function isNotBelow{T<:Real}(value::T, bounds::BotTop{T})
    !isBelow(value,upperBound)
end


      
      

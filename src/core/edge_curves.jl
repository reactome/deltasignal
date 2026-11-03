# specs/128: response curves on edges (specs/125 stage 1).
#
# Each reaction slot reads its input's fold f after every structural rule has decided what it reads; the curve decides
# how strongly the reaction responds. Activating slots: g(f) = f^h (1 + K^h) / (f^h + K^h), so g(1) = 1, g(0) = 0, the
# elasticity at rest is h K^h / (1 + K^h) and the most the edge passes up is 1 + K^h; K = Inf, h = 1 is today's
# identity. Inhibiting slots (negative regulator under `divide`, depletion): (1 + k^m) / (1 + k^m f^m), the 2025
# design's 1 / (1 + (I/K_I)^m) normalised to 1 at rest. Output edges and identity links are never curved.

# Slot classes, set from the edge type in `create_reaction_from_edges`.
const SLOT_IDENTITY = 0x00
const SLOT_INPUT = 0x01        # input, assembly, composition
const SLOT_CATALYST = 0x02
const SLOT_REGULATOR_POS = 0x03

const CURVE_ACTIVATING = ("input", "catalyst", "regulator_pos")
const CURVE_INHIBITING = ("regulator_neg", "depletion")

function slot_class(edge_type::AbstractString)::UInt8
    edge_type in ("input", "assembly", "composition") && return SLOT_INPUT
    edge_type == "catalyst" && return SLOT_CATALYST
    edge_type == "regulator" && return SLOT_REGULATOR_POS
    return SLOT_IDENTITY
end

"""
The edge-curve table in force. `K`/`h` are indexed by slot class (input, catalyst, regulator_pos); an activating class
at K = Inf and h = 1 is the identity and is skipped. `inh`/`dep` are `nothing` for today's divide form, or (k, m).
"""
struct EdgeCurves
    on::Bool
    K::NTuple{3, Float64}
    h::NTuple{3, Float64}
    inh::Union{Nothing, Tuple{Float64, Float64}}
    dep::Union{Nothing, Tuple{Float64, Float64}}
    desc::String
end

const EDGE_CURVES_OFF = EdgeCurves(false, (Inf, Inf, Inf), (1.0, 1.0, 1.0), nothing, nothing, "off")
const _EDGE_CURVE_CACHE = Dict{Tuple{String, Float64}, EdgeCurves}()

curve_active(c::EdgeCurves, cls::UInt8) = c.on && cls != SLOT_IDENTITY && !(isinf(c.K[cls]) && c.h[cls] == 1.0)

"""Activating curve on a fold (f >= 0). Identity when K = Inf and h = 1."""
function curve_fold(f::T, K::Float64, h::Float64)::T where {T<:Real}
    f <= zero(T) && return zero(T)
    fh = h == 1.0 ? f : f^T(h)
    isinf(K) && return fh
    Kh = K^h
    return fh * T(1 + Kh) / (fh + T(Kh))
end

"""Inhibiting factor for an inhibitor at fold f (f >= 0): (1 + k^m) / (1 + k^m f^m)."""
function curve_inhibit(f::T, k::Float64, m::Float64)::T where {T<:Real}
    km = k^m
    fm = f <= zero(T) ? zero(T) : (m == 1.0 ? f : f^T(m))
    return T(1 + km) / (one(T) + T(km) * fm)
end

function _curve_number(raw::AbstractString, what::String; inf_ok::Bool=false)::Float64
    s = lowercase(strip(raw))
    inf_ok && s == "inf" && return Inf
    v = tryparse(Float64, s)
    (v === nothing || !(v > 0) || !isfinite(v)) &&
        throw(ArgumentError("DS_EDGE_CURVES: $what = \"$raw\" must be a positive number" * (inf_ok ? " or inf" : "")))
    return v
end

"""
Parse the table at `path`: a header `class	K	h`, then one row per class. Activating classes take K (a number or
`inf`) and h; inhibiting classes take k (a number, or `divide` for today's form) and m. A class not listed keeps
today's behaviour.
"""
function parse_edge_curves(path::AbstractString)::EdgeCurves
    isfile(path) || throw(ArgumentError("DS_EDGE_CURVES=$path is not a file"))
    lines = [l for l in readlines(path) if !isempty(strip(l)) && !startswith(strip(l), "#")]
    isempty(lines) && throw(ArgumentError("DS_EDGE_CURVES=$path is empty"))
    split(lines[1], '\t') == ["class", "K", "h"] ||
        throw(ArgumentError("DS_EDGE_CURVES: the header must be class<TAB>K<TAB>h, got \"$(lines[1])\""))
    K = [Inf, Inf, Inf]; h = [1.0, 1.0, 1.0]
    inh = nothing; dep = nothing
    seen = Set{String}()
    for l in lines[2:end]
        f = split(l, '\t')
        length(f) == 3 || throw(ArgumentError("DS_EDGE_CURVES: \"$l\" must have three tab-separated fields"))
        cls = strip(f[1])
        cls in seen && throw(ArgumentError("DS_EDGE_CURVES: class $cls is listed twice"))
        push!(seen, cls)
        if cls in CURVE_ACTIVATING
            i = findfirst(==(cls), CURVE_ACTIVATING)
            K[i] = _curve_number(f[2], "$cls K"; inf_ok=true)
            h[i] = _curve_number(f[3], "$cls h")
        elseif cls in CURVE_INHIBITING
            km = if lowercase(strip(f[2])) == "divide"
                lowercase(strip(f[3])) in ("divide", "1") ||
                    throw(ArgumentError("DS_EDGE_CURVES: $cls is divide, so its m must be 1 or divide"))
                nothing
            else
                (_curve_number(f[2], "$cls k"), _curve_number(f[3], "$cls m"))
            end
            cls == "regulator_neg" ? (inh = km) : (dep = km)
        else
            throw(ArgumentError("DS_EDGE_CURVES: unknown class \"$cls\"; classes are " *
                                join((CURVE_ACTIVATING..., CURVE_INHIBITING...), ", ")))
        end
    end
    desc = join(["$(c) K=$(K[i]) h=$(h[i])" for (i, c) in enumerate(CURVE_ACTIVATING)], "; ") *
           "; regulator_neg " * (inh === nothing ? "divide" : "k=$(inh[1]) m=$(inh[2])") *
           "; depletion " * (dep === nothing ? "divide" : "k=$(dep[1]) m=$(dep[2])")
    return EdgeCurves(true, Tuple(K), Tuple(h), inh, dep, desc)
end

"""`DS_EDGE_CURVES`: `off` (default) or the path of a curve table; parsed once per file version."""
function edge_curves_env()::EdgeCurves
    raw = strip(get(ENV, "DS_EDGE_CURVES", ""))
    (isempty(raw) || lowercase(raw) == "off") && return EDGE_CURVES_OFF
    key = (String(raw), isfile(raw) ? mtime(raw) : 0.0)
    return get!(() -> parse_edge_curves(raw), _EDGE_CURVE_CACHE, key)
end

"""specs/128 report: the number of network edges in each class whose curve is not the identity."""
function edge_curve_slots(network, c::EdgeCurves)::Dict{String, Int}
    out = Dict{String, Int}()
    for e in network.edges
        name = if e.edge_type == "depletion"
            c.dep === nothing ? nothing : "depletion"
        elseif !e.is_positive
            c.inh === nothing ? nothing : "regulator_neg"
        else
            cls = slot_class(e.edge_type)
            curve_active(c, cls) ? CURVE_ACTIVATING[cls] : nothing
        end
        name === nothing || (out[name] = get(out, name, 0) + 1)
    end
    return out
end

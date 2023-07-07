using CairoMakie
using CoherentNoise
using Images

##

struct BiomeRegion
    tempMin::Float64
    tempMax::Float64
    humidityMin::Float64
    humidityMax::Float64
    biome::Symbol
end

biomecolors = Dict(
    :desert => colorant"yellow",
    :ocean => colorant"blue",
    :jungle => colorant"darkgreen",
    :ice => colorant"lightcyan",
    :mountain => colorant"lightgrey",
    :savanna => colorant"darkolivegreen",
    :plains => colorant"lightgreen",
    :forest => colorant"green",
    :beach => colorant"tan",
    :unknown => colorant"grey",
    :pink => colorant"deeppink",
    :purple => colorant"purple"
)

function distance(biome::BiomeRegion, x, y)
    if (x >= biome.tempMin && x <= biome.tempMax)
        return min(abs(y - biome.humidityMax), abs(y - biome.humidityMin))
    end
    if (y >= biome.humidityMin && y <= biome.humidityMax)
        return min(abs(x - biome.tempMax), abs(x - biome.tempMin))
    end
    return min(
        euclidean([x, y], (biome.tempMax, biome.humidityMax)),
        euclidean([x, y], (biome.tempMax, biome.humidityMin)),
        euclidean([x, y], (biome.tempMin, biome.humidityMax)),
        euclidean([x, y], (biome.tempMin, biome.humidityMin))
    )
end

function paramstobiome(x, y, biomes::Vector{BiomeRegion}; extrapolate = true)
    min = sqrt(2)
    out = :unknown
    for biome in biomes
        if biome.tempMin <= x <= biome.tempMax && biome.humidityMin <= y <= biome.humidityMax
            out = biome.biome
            break
        end
        if extrapolate
            dist = distance(biome, x, y)
            if dist < min
                min = dist
                out = biome.biome
            end
        end
    end
    out
end

function visualizebiomes(biomes::Vector{BiomeRegion}; resolution = 512, extrapolate = true)
    img = zeros(RGB{Float64}, resolution, resolution)
    for i in 1:resolution
        for j in 1:resolution
            x = i / resolution
            y = j / resolution
            img[i, j] = biomecolors[paramstobiome(x, y, biomes, extrapolate = extrapolate)]
        end
    end
    img
end

function visualizebiomes(layers::Vector{Vector{BiomeRegion}}; resolution = 512, extrapolate = true)
    img = zeros(RGB{Float64}, resolution, resolution)
    for i in 1:resolution
        for j in 1:resolution
            x = i / resolution
            y = j / resolution
            biome = :unknown
            for layer in layers
                biome = paramstobiome(x, y, layer, extrapolate = false)
                if biome != :unknown
                    break
                end
            end
            if extrapolate && biome == :unknown
                biome = paramstobiome(x, y, layers[end], extrapolate = true)
            end
            img[i, j] = biomecolors[biome]
        end
    end
    img
end

##

originalbiomes = [
    BiomeRegion(0.1, 1.0, 0.6, 1.0, :ocean)
    BiomeRegion(0.7, 1.0, 0.0, 0.4, :desert)
    BiomeRegion(0.6, 1.0, 0.4, 0.55, :jungle)
    BiomeRegion(0.0, 0.1, 0.5, 1.0, :ice)
    BiomeRegion(0.0, 0.1, 0.0, 0.5, :mountain)
    BiomeRegion(0.1, 0.7, 0.0, 0.2, :savanna)
    BiomeRegion(0.1, 0.7, 0.2, 0.4, :plains)
    BiomeRegion(0.1, 0.6, 0.4, 0.55, :forest)
    BiomeRegion(0.1, 1.0, 0.55, 0.6, :beach)
]

keylist = [b.biome for b in originalbiomes]

function biomeregionsgenerate(size, sampler; seed = 42, scale = 1)
    sampler1 = opensimplex2_2d(seed = seed^7 % 5) * 0.5 + 0.5
    sampler2 = opensimplex2_2d(seed = seed^11 % 3) * 0.5 + 0.5
    out = fill(:nothing, size, size)
    for i ∈ 1:size
        for j ∈ 1:size
            x = i * scale * 0.005
            y = j * scale * 0.005
            out[i, j] = sampler(sample(sampler1, x, y), sample(sampler2, x, y))
        end
    end
    out
end

function regionsampler(biomes)
    (x, y) -> paramstobiome(x, y, biomes)
end

function biomestoimg(biomemap)
    img = zeros(RGB{Float64}, size(biomemap))
    for i in eachindex(biomemap)
        img[i] = biomecolors[biomemap[i]]
    end
    img
end

function generatedbiomeratios(biomemap)
    out = Dict{Symbol, Float64}()
    count = 0
    for biome ∈ biomemap
        if !haskey(out, biome)
            out[biome] = 0.0
        end
        out[biome] += 1
        count += 1
    end
    for k in keys(out)
        out[k] /= count
    end
    out
end

function withreplacementsampler(layers)
    (x, y) -> begin
        for layer in layers
            sample = paramstobiome(x, y, layer, extrapolate = false)
            if sample != :unknown
                return sample
            end
        end
        return paramstobiome(x, y, layers[end], extrapolate = true)
    end
end

overlaybiomes = [
    BiomeRegion(0.2, 0.4, 0.1, 0.3, :pink)
    BiomeRegion(0.8, 1.0, 0.5, 0.7, :purple)
]

##

if !isdir("out")
    mkdir("out")
end

visualizeoriginal = visualizebiomes(originalbiomes)
exampleoriginal = biomestoimg(biomeregionsgenerate(512, regionsampler(originalbiomes)))

save("out/visualizebiomes.png", visualizeoriginal)
save("out/exampleoriginal.png", exampleoriginal)

visualizeoverlay = visualizebiomes([overlaybiomes, originalbiomes])
exampleoverlay = biomestoimg(biomeregionsgenerate(512, withreplacementsampler([overlaybiomes, originalbiomes])))

save("out/visualizeoverlay.png", visualizeoverlay)
save("out/exampleoverlay.png", exampleoverlay)

##

originalrarities = generatedbiomeratios(biomeregionsgenerate(4096, regionsampler(originalbiomes)))

##

f = Figure()
ax = Makie.Axis(
    f[1,1],
    xticks = (1:length(keylist), string.(keylist)),
    ylabel = "Biome Rarity Relative to Original"
)
newrarities = generatedbiomeratios(biomeregionsgenerate(4096, withreplacementsampler([overlaybiomes, originalbiomes])))

biomescarcitychange = Dict()
for k in keys(originalrarities)
    biomescarcitychange[k] = newrarities[k] / originalrarities[k]
end

for i in 1:length(keylist)
    plot!(
        ax,
        [i],
        [biomescarcitychange[keylist[i]]],
        color = biomecolors[keylist[i]],
        markersize = 50,
        strokecolor = :black,
        strokewidth = 5
    )
end

save("out/overlayrarities.png", f)
f

##

function visualizebiomesampler(sampler; resolution = 512)
    img = zeros(RGB{Float64}, resolution, resolution)
    for i in 1:resolution
        for j in 1:resolution
            x = i / resolution
            y = j / resolution
            img[i, j] = biomecolors[sampler(x, y)]
        end
    end
    img
end

function distancetoedge(x, y)
    min(0.5 - abs(x - 0.5), 0.5 - abs(y - 0.5))
end

function distancetoedgeprojection(x, y, point)
    if (x == point[1] && y == point[2])
        return distancetoedge(x, y)
    end
    if (x <= 0 || y <= 0 || x >= 1.0 || y >= 1.0)
        return 0
    end
    topoint = point .- [x, y]
    if topoint[1] == 0
        return if topoint[2] < 0 1 - y else y end
    end
    if topoint[2] == 0
        return if topoint[1] < 0 1 - x else x end
    end

    yInter1 = y - x / topoint[1] * topoint[2]
    xInter1 = x - y / topoint[2] * topoint[1]
    yInter2 = y + (1 - x) / topoint[1] * topoint[2]
    xInter2 = x + (1 - y) / topoint[2] * topoint[1]

    if (topoint[2] > 0 && xInter1 >= 0 && xInter1 <= 1)
        euclidean((xInter1, 0), (x, y))
    elseif (topoint[2] < 0 && xInter2 >= 0 && xInter2 <= 1)
        euclidean((xInter2, 1), (x, y))
    elseif (topoint[1] > 0 && yInter1 >= 0 && yInter1 <= 1)
        euclidean((0, yInter1), (x, y))
    elseif (topoint[1] < 0 && yInter2 >= 0 && yInter2 <= 1)
        euclidean((1, yInter2), (x, y))
    else
        0
    end
end

function transformedsampler(sampler, point, cutoff)
    (x, y) -> begin
        edgepointdist = distancetoedgeprojection(x, y, point)
        pointdist = euclidean((x, y), point)
        dist = pointdist/(edgepointdist+pointdist)
        if dist < cutoff
            return :unknown
        end
        movedistratio = (dist - sqrt((dist^2 - cutoff^2)/(1 - cutoff^2))) / cutoff
        topoint = point .- [x, y]
        topoint .*= movedistratio
        sampler(x + topoint[1], y + topoint[2])
    end
end

##

function addregion(sampler, region, cumulativearea)
    center = [(region.tempMax + region.tempMin)/2, (region.humidityMax + region.humidityMin)/2]
    area = (region.tempMax - region.tempMin) * (region.humidityMax - region.humidityMin)
    area /= (cumulativearea + area)
    sidelength = sqrt(area)
    sampler = transformedsampler(sampler, center, sidelength)
    lowerleft = center .* (1 - sidelength)
    bottomright = (center .- 1) .* (1 - sidelength) .+ 1
    (x, y) -> begin
        if x < lowerleft[1] || y < lowerleft[2] || x > bottomright[1] || y > bottomright[2]
            return sampler(x, y)
        else
            return region.biome
        end
    end, cumulativearea + area
end

##

sampler = regionsampler(originalbiomes)
cumulativearea = 1
for region in overlaybiomes
    sampler, cumulativearea = addregion(sampler, region, cumulativearea)
end

visualizestretched = visualizebiomesampler(sampler)
examplestretched = biomestoimg(biomeregionsgenerate(512, sampler))

save("out/visualizestretched.png", visualizestretched)
save("out/examplestretched.png", examplestretched)

##

f = Figure()
ax = Makie.Axis(
    f[1,1],
    xticks = (1:length(keylist), string.(keylist)),
    ylabel = "Biome Rarity Relative to Original"
)
newrarities = generatedbiomeratios(biomeregionsgenerate(4096, sampler))

biomescarcitychange = Dict()
for k in keys(originalrarities)
    biomescarcitychange[k] = newrarities[k] / originalrarities[k]
end

for i in 1:length(keylist)
    plot!(
        ax,
        [i],
        [biomescarcitychange[keylist[i]]],
        color = biomecolors[keylist[i]],
        markersize = 50,
        strokecolor = :black,
        strokewidth = 5
    )
end

save("out/stretchedrarities.png", f)
f

##
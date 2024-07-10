
function produce_video(dir)
   imgnames = filter(x->occursin(".png",x), readdir(dir)) # Populate list of all .pngs
   namestrings =  map(x->split(x,".")[1], imgnames) # Extract index from filenames
   intstrings =  map(x->split(x,"=")[2], namestrings) # Extract index from filenames
   p = sortperm(parse.(Int, intstrings)) #sort files numerically
   imgnames = imgnames[p]

   encoder_options = (crf=23, preset="medium")


   imgstack=Array{Matrix{RGB{N0f8}}}(undef, length(imgnames))

   for i in collect(1:length(imgnames))
      imgstack[i]=load(joinpath(dir, imgnames[i]))
   end


   VideoIO.save(dir*"video.mp4", imgstack, framerate=1, encoder_options=encoder_options)
end


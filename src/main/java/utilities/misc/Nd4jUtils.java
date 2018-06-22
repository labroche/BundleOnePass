package utilities.misc;

import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.api.ops.impl.accum.Max;
import org.nd4j.linalg.api.ops.impl.accum.Min;
import org.nd4j.linalg.cpu.nativecpu.NDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import java.util.stream.Collectors;

public class Nd4jUtils {
    private static final Logger LOGGER = Logger.getLogger( Nd4jUtils.class.getName() );

    public static INDArray normalize(INDArray data) {
        INDArray out = Nd4j.create(data.shape());
        for (int dim = 0; dim < out.columns(); dim++) {
            double min = Nd4j.getExecutioner().execAndReturn(new Min(data.getColumn(dim))).getFinalResult().doubleValue();
            double max = Nd4j.getExecutioner().execAndReturn(new Max(data.getColumn(dim))).getFinalResult().doubleValue();
            double delta = max - min;
            out.putColumn(dim, data.getColumn(dim).sub(min).div(delta));
        }
        return out;
    }

    public static void normalizeInPlace(INDArray data) {
        for (int dim = 0; dim < data.columns(); dim++) {
            double min = Nd4j.getExecutioner().execAndReturn(new Min(data.getColumn(dim))).getFinalResult().doubleValue();
            double max = Nd4j.getExecutioner().execAndReturn(new Max(data.getColumn(dim))).getFinalResult().doubleValue();
            double delta = max - min;
            data.getColumn(dim).subi(min).divi(delta);
        }
    }

    public static ArrayList<String> loadIds(String path, String separator, int pos){
        try {
            return Files.lines(Paths.get(path)).map(line -> line.split(separator)[pos]).collect(Collectors.toCollection(ArrayList::new));
        } catch (IOException e) {
            LOGGER.warning("Cannot find file at '" + path + "'.");
            return new ArrayList<>();
        }
    }

    public static ArrayList<String> loadIds(String path, String separator, int pos, String toDelete){
        try {
            return Files.lines(Paths.get(path)).map(line -> line.split(separator)[pos].replace(toDelete, "")).collect(Collectors.toCollection(ArrayList::new));
        } catch (IOException e) {
            LOGGER.warning("Cannot find file at '" + path + "'.");
            return new ArrayList<>();
        }
    }

    public static INDArray loadCSV(String path, String separator, int label, boolean ignoreFirstLine){
        INDArray out;
        ArrayList<Float[]> temp = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;
            if(ignoreFirstLine)
                br.readLine();
            while ((line = br.readLine()) != null) {
                String[] cols = line.split(separator);
                ArrayList<Float> thisline = new ArrayList<>();
                for(int i = 0; i < cols.length; ++i){
                    if (i != label){
                        try {
                            thisline.add(Float.parseFloat(cols[i]));
                        }catch (NumberFormatException e){
                            System.out.printf("ERROR Line '%s'%n", thisline);
                        }

                    }
                }
                temp.add(thisline.toArray(new Float[0]));
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        float[][] tmp = new float[temp.size()][temp.get(0).length];
        for (int i = 0; i < temp.size(); ++i){
            for (int j = 0; j < temp.get(0).length; ++j)
                tmp[i][j] = temp.get(i)[j];
        }
        out = new NDArray(tmp);
        //System.out.println(out);
        return out;
    }

    public static List<INDArray> loadCSVAsList(String path, String separator, int label, boolean ignoreFirstLine){
        List<INDArray> out = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;
            if(ignoreFirstLine)
                br.readLine();
            while ((line = br.readLine()) != null) {
                String[] cols = line.split(separator);
                float[] tmp = new float[label == -1 ? cols.length : cols.length - 1];
                for (int i = 0; i < tmp.length; i++) {
                    if (i != label)
                        tmp[i] = Float.parseFloat(cols[i]);
                }
                out.add(Nd4j.create(tmp));
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        return out;
    }

    public static List<float[]> loadCSVAsListNative(String path, String separator, int label, boolean ignoreFirstLine){
        List<float[]> out = new ArrayList<>();
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;
            if(ignoreFirstLine)
                br.readLine();
            while ((line = br.readLine()) != null) {
                String[] cols = line.split(separator);
                float[] tmp = new float[label == -1 ? cols.length : cols.length - 1];
                for (int i = 0; i < tmp.length; i++) {
                    if (i != label)
                        tmp[i] = Float.parseFloat(cols[i]);
                }
                out.add(tmp);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return out;
    }

    public static INDArray toMatrix(List<INDArray> input){
        INDArray out = Nd4j.create(input.size(), input.get(0).columns());
        int i = 0;
        for (INDArray x : input) {
            out.putRow(i++, x);
        }
        return out;
    }

    public static <E> void writePartition(List<List<E>> partition, String path) throws FileNotFoundException{
        File f = new File(path);
        if (f.exists())
            f.delete();
        PrintWriter out = new PrintWriter(new FileOutputStream(f));
        int i = 0;
        for (List<E> cluster : partition){
            for (E element : cluster){
                out.println(element.toString() + ";" + i);
            }
            i++;
        }
        out.close();
    }
}

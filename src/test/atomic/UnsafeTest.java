/**
 * @program: LearningJDK
 * @description:
 * @author: Shaobo.Qian
 * @create: 2019-12-19 15:04
 **/

package test.atomic;

import sun.misc.Unsafe;

public class UnsafeTest {
    public static void main(String[] args) {
        Unsafe unsafe = Unsafe.getUnsafe();
        System.out.println("unsafe = " + unsafe);
    }

}

"""核对配对WLS、HC3、抵消均值及错误率分母的基础测试。"""
import unittest
import numpy as np
from pipeline.calibration import geometry,evaluate,mean_profiles,wilson


class PairedTests(unittest.TestCase):
    """用独立矩阵表达核查分析核，不以预期显著性作为通过条件。"""
    def setUp(self):
        """每个测试使用同一可复现、非退化的嵌套设计。"""
        self.rng=np.random.default_rng(41)
        self.x0=np.column_stack((np.ones(60),self.rng.integers(0,2,60),self.rng.normal(size=60)))
        self.x1=np.column_stack((self.x0,self.rng.normal(size=(60,3))))
        self.w=np.exp(self.rng.normal(0,.3,60))
        self.g=geometry(self.x0,self.x1,self.w)

    def test_coefficients(self):
        """正规方程映射应与加权SVD最小二乘一致。"""
        y=self.rng.normal(size=(60,4));d,*_=evaluate(y,self.g)
        b0=np.linalg.lstsq(self.x0*np.sqrt(self.w[:,None]),y*np.sqrt(self.w[:,None]),rcond=None)[0]
        b1=np.linalg.lstsq(self.x1*np.sqrt(self.w[:,None]),y*np.sqrt(self.w[:,None]),rcond=None)[0]
        np.testing.assert_allclose(d,b1[1]-b0[1],atol=1e-12)

    def test_covariance(self):
        """同响应差值方差应等于完整模型中对应线性对比的方差。"""
        g=self.g;c=np.eye(self.x1.shape[1])[1]-g['a0'][1]@self.x1
        v=c@np.linalg.inv(self.x1.T@(self.w[:,None]*self.x1))@c
        self.assertAlmostEqual(v,g['unit'],12)

    def test_hc3(self):
        """直接sandwich矩阵与观测加权求和两种HC3表达应一致。"""
        g=self.g;y=self.rng.normal(size=(60,1));_,_,se,_=evaluate(y,g)
        e=(y-self.x1@(g['a1']@y)).ravel()
        cov=g['a1']@np.diag(e*e/(1-g['h'])**2)@g['a1'].T
        c=np.eye(self.x1.shape[1])[1]-g['a0'][1]@self.x1
        self.assertAlmostEqual(se[0]**2,c@cov@c,12)

    def test_means(self):
        """质量项可非零而目标差值为零，避免把F检验误当目标变化检验。"""
        means=mean_profiles(self.g)
        self.assertAlmostEqual(self.g['d']@means[0][1],0,12)
        self.assertAlmostEqual(self.g['d']@means[1][1],0,12)
        self.assertGreater(np.linalg.norm(means[1][1]-means[0][1]),.1)
        self.assertAlmostEqual(self.g['d']@means[2][1],2*np.sqrt(self.g['unit']),12)

    def test_degenerate(self):
        """质量列与基础空间正交时系数差恒为零，不能构造常规差异检验。"""
        x0=self.x0;z=self.x1[:,-3:];z=z-x0@np.linalg.lstsq(x0,z,rcond=None)[0]
        g=geometry(x0,np.column_stack((x0,z)),np.ones(60))
        self.assertLess(g['unit'],1e-14)
        with self.assertRaises(ValueError):mean_profiles(g)

    def test_wilson(self):
        """零事件和全事件仍有非零宽度的蒙特卡洛区间。"""
        self.assertGreater(wilson(0,2500)[1],0)
        self.assertLess(wilson(2500,2500)[0],1)

    def test_non_nested(self):
        """非嵌套矩阵不得静默进入本检验。"""
        with self.assertRaises(ValueError):geometry(self.x0,self.rng.normal(size=(60,6)),self.w)


if __name__=='__main__':unittest.main(verbosity=2)
